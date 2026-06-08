/*!
 * \file floe/lcp/solver/GS_solver.hpp
 * \brief Projected Gauss-Seidel (NSCD-style) contact solver — alternative to the Lemke solver.
 * \author Quentin Jouet
 *
 * This is a *second* contact solver, kept fully separate from the historical Lemke pivoting solver
 * (LCP_solver.hpp). It is meant to be routed *only* on large, anchored, quasi-static contact
 * components (dense jams), where the Lemke + active-subgraph propagation explodes in number of LCP
 * and fails to converge (hyperstatic / degenerate systems). On such components a single projected
 * Gauss-Seidel solve of the *whole* component:
 *   - converges where pivoting cycles (degeneracy-tolerant: it relaxes to a feasible force network
 *     instead of pivoting to a unique vertex),
 *   - has a deterministic O(contacts) per-sweep cost (no pivot-count spikes),
 *   - needs no random perturbation of the system,
 *   - leaves no active contact afterwards, so the adaptive time step can stay large.
 *
 * It reuses the existing assembly (GraphLCP) for the mass matrix, the inverse mass matrix, the
 * normal/tangential Jacobians and the friction coefficients, and reuses GraphLCP::apply_impulses to
 * feed the contact impulses back into the graph exactly like the Lemke path. The contact law solved
 * here is Signorini (non-penetration) + Coulomb friction, purely inelastic (restitution e = 0),
 * which is the relevant law for held/jammed contacts. Restitution (for genuine impulsive collisions,
 * e.g. Newton's cradle) is intentionally out of scope here: those stay on the Lemke solver.
 */

#ifndef FLOE_LCP_SOLVER_GS_SOLVER_HPP
#define FLOE_LCP_SOLVER_GS_SOLVER_HPP

#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/range/iterator_range.hpp>

#include "floe/lcp/builder/graph_to_lcp.hpp"
#include "../product/config/config_base.hpp" // types

namespace floe { namespace lcp { namespace solver
{

namespace ublas = boost::numeric::ublas;

/*! Projected Gauss-Seidel contact solver (one solve over a whole contact component).
 *
 * \tparam T Fundamental scalar type.
 */
template <typename T>
class GaussSeidelSolver
{
public:
    GaussSeidelSolver(int max_iter = 2000, T tol = T(1e-8))
        : m_max_iter(max_iter), m_tol(tol) {}

    void set_params(int max_iter, T tol) {
        if (max_iter > 0) m_max_iter = max_iter;
        if (tol > 0)      m_tol = tol;
    }
    void set_max_iter(int max_iter) { if (max_iter > 0) m_max_iter = max_iter; }

    int  max_iter() const { return m_max_iter; }
    T    tol()      const { return m_tol; }

    /*! Solve the contact problem on the whole component \p graph by projected Gauss-Seidel.
     *
     *  Writes back the post-contact velocities into the floes and the contact impulses into the graph
     *  (via GraphLCP::apply_impulses), exactly like the Lemke path.
     *
     *  \param[in]  graph         A connected contact component (sub-graph).
     *  \param[out] out_iters     Number of sweeps performed.
     *  \param[out] out_residual  Final residual (max penetration velocity reached).
     *  \param[out] out_vmax_free Max translational free (driving) speed of the component — used by the
     *                            caller to decide, per floe, which floes the solve actually brought to rest.
     *  \return true if a solution was committed (best-effort); false if it declined to genuine garbage.
     */
    template <typename TGraph>
    bool solve(TGraph const& graph, int& out_iters, T& out_residual, T& out_vmax_free) const
    {
        using floe::lcp::builder::GraphLCP;

        GraphLCP<T, TGraph> glcp(graph); // assembles M, invM, J, D, mu (no big LCP matrix built)
        const int n = glcp.nb_floes;
        const int m = glcp.nb_contacts;
        out_iters = 0;
        out_residual = 0;
        out_vmax_free = 0;
        if (m == 0) return true;

        // Free velocity v = current floe velocities (already include the wind/ocean drag of the
        // previous move; solve_contacts runs before the next move). This is the "restart" velocity.
        ublas::vector<T> v(3 * n);
        for (int i = 0; i < n; ++i) {
            auto const st = graph[i].floe->state();
            v(3 * i)     = st.speed.x;
            v(3 * i + 1) = st.speed.y;
            v(3 * i + 2) = st.rot;
        }

        // Per-contact local data: the (up to) 6 dof rows (3 per floe), the normal and tangential
        // jacobian values, the inverse-mass on each row, and the diagonal Delassus blocks Wnn / Wtt.
        struct Row { int idx; T jn; T dtg; T im; };
        std::vector<std::array<Row, 6>> rows(m);
        std::vector<T> Wnn(m, 0), Wtt(m, 0), muv(m, 0);
        {
            int a = 0;
            for (auto const& edge : boost::make_iterator_range(edges(graph)))
            {
                const int s  = static_cast<int>(source(edge, graph));
                const int t  = static_cast<int>(target(edge, graph));
                const int b1 = 3 * std::min(s, t);
                const int b2 = 3 * std::max(s, t);
                const int base[6] = { b1, b1 + 1, b1 + 2, b2, b2 + 1, b2 + 2 };
                for (std::size_t c = 0; c < graph[edge].size(); ++c)
                {
                    for (int k = 0; k < 6; ++k)
                    {
                        const int idx = base[k];
                        const T jn  = glcp.J(idx, a);
                        const T dtg = glcp.D(idx, 2 * a);
                        const T im  = glcp.invM(idx, idx);
                        rows[a][k] = Row{ idx, jn, dtg, im };
                        Wnn[a] += im * jn * jn;
                        Wtt[a] += im * dtg * dtg;
                    }
                    muv[a] = glcp.mu(a, a);
                    ++a;
                }
            }
        }

        // Contact impulses (cold start; warm-start across steps is a planned next step).
        ublas::vector<T> pn(m, 0), pt(m, 0);

        // Relative scale used to harden the degenerate-block guard against division blow-ups, and to
        // make the convergence test and the divergence cap *scale-relative* (FloeDyn runs span domains
        // scaled by orders of magnitude, so absolute thresholds are meaningless).
        T Wnn_max = 0, vmax_free = 0;
        for (int a = 0; a < m; ++a) Wnn_max = std::max(Wnn_max, Wnn[a]);
        for (int i = 0; i < n; ++i) // translational free speed of each floe
            vmax_free = std::max(vmax_free, std::sqrt(v(3*i)*v(3*i) + v(3*i+1)*v(3*i+1)));
        out_vmax_free = vmax_free;
        const T W_guard = std::max(T(1e-30), Wnn_max * T(1e-12));
        const T tol_feasible = std::max(T(1e-9), m_tol * vmax_free); // "converged" reporting threshold

        // Keep the BEST (lowest penetration velocity) solution seen, and commit it even if not perfectly
        // converged. FloeDyn philosophy: a slightly imperfect contact solution that lets the simulation
        // advance beats an exact one we can never reach (especially with --bustle 0, where Lemke struggles
        // too). We only decline (and fall back to Lemke) on genuine numerical garbage.
        ublas::vector<T> v_best = v, pn_best = pn, pt_best = pt;
        T pen_best = std::numeric_limits<T>::max();
        bool finite_ok = true;

        int it = 0;
        bool converged = false;
        for (; it < m_max_iter; ++it)
        {
            T pen = 0;
            for (int a = 0; a < m; ++a)
            {
                if (Wnn[a] <= W_guard) continue; // contact between (near-)fixed bodies

                // --- Normal: Signorini, inelastic (target normal velocity 0) ---
                T un = 0;
                for (int k = 0; k < 6; ++k) un += rows[a][k].jn * v(rows[a][k].idx);
                if (!std::isfinite(un)) { finite_ok = false; break; }
                if (-un > pen) pen = -un; // feasibility = no approaching (penetrating) relative velocity

                T pn_new = pn(a) - un / Wnn[a];
                if (pn_new < 0) pn_new = 0;
                if (!std::isfinite(pn_new)) { finite_ok = false; break; }
                const T dpn = pn_new - pn(a);
                if (dpn != T(0)) {
                    for (int k = 0; k < 6; ++k) { auto const& r = rows[a][k]; v(r.idx) += r.im * r.jn * dpn; }
                    pn(a) = pn_new;
                }

                // --- Tangential: Coulomb friction, |pt| <= mu * pn ---
                if (Wtt[a] > W_guard)
                {
                    T ut = 0;
                    for (int k = 0; k < 6; ++k) ut += rows[a][k].dtg * v(rows[a][k].idx);
                    if (!std::isfinite(ut)) { finite_ok = false; break; }
                    T pt_new = pt(a) - ut / Wtt[a];
                    const T bound = muv[a] * pn(a);
                    if (pt_new >  bound) pt_new =  bound;
                    if (pt_new < -bound) pt_new = -bound;
                    const T dpt = pt_new - pt(a);
                    if (dpt != T(0)) {
                        for (int k = 0; k < 6; ++k) { auto const& r = rows[a][k]; v(r.idx) += r.im * r.dtg * dpt; }
                        pt(a) = pt_new;
                    }
                }
            }
            if (!finite_ok) break; // diverged to non-finite: stop and keep the best earlier solution

            if (pen < pen_best) { pen_best = pen; v_best = v; pn_best = pn; pt_best = pt; }
            if (pen < tol_feasible) { converged = true; ++it; break; }
        }
        out_iters = it;
        out_residual = pen_best;

        // Decline (-> Lemke fallback) only on genuine garbage: no usable solution at all, or a best
        // solution that is non-finite or diverged far beyond any physical velocity.
        if (pen_best == std::numeric_limits<T>::max()) return false;
        T vmax_post = 0;
        for (int i = 0; i < n; ++i)
            vmax_post = std::max(vmax_post, std::sqrt(v_best(3*i)*v_best(3*i) + v_best(3*i+1)*v_best(3*i+1)));
        const T V_CAP = std::max(T(1e3), T(100) * vmax_free);
        if (!std::isfinite(vmax_post) || vmax_post > V_CAP) return false;
        (void)converged;

        // Commit the best velocities (obstacles keep their state).
        for (int i = 0; i < n; ++i)
        {
            if (graph[i].floe->is_obstacle()) continue;
            auto& st = graph[i].floe->state();
            st.speed = { v_best(3 * i), v_best(3 * i + 1) };
            st.rot   = v_best(3 * i + 2);
        }

        // Feed the best contact impulses back into the graph (same bookkeeping as the Lemke path).
        ublas::vector<T> normal(m), tangential(2 * m);
        for (int a = 0; a < m; ++a) {
            normal(a)          = pn_best(a);
            tangential(2 * a)     = pt_best(a);
            tangential(2 * a + 1) = 0;
        }
        glcp.apply_impulses(normal, tangential);

        return true; // committed (best-effort)
    }

private:
    int m_max_iter; //!< maximum number of sweeps
    T   m_tol;      //!< convergence tolerance (velocity-scaled impulse change of a sweep)
};

}}} // namespace floe::lcp::solver

#endif // FLOE_LCP_SOLVER_GS_SOLVER_HPP
