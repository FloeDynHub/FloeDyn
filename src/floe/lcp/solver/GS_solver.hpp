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
#include <map>
#include <utility>

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
    GaussSeidelSolver(int max_iter = 20000, T tol = T(1e-8))
        : m_max_iter(max_iter), m_tol(tol) {}

    void set_params(int max_iter, T tol) {
        if (max_iter > 0) m_max_iter = max_iter;
        if (tol > 0)      m_tol = tol;
    }
    void set_max_iter(int max_iter) { if (max_iter > 0) m_max_iter = max_iter; }

    int  max_iter() const { return m_max_iter; }
    T    tol()      const { return m_tol; }

    /*! Two solve modes (the "two passes" of OPTIMJAM, run per component per step):
     *  - Forces:   true masses for every floe (only real obstacles are fixed), reads the floes' free
     *              (drag-loaded) pre-resolution velocities — which is what sets the pressure on the
     *              boundary. RECORDS the resulting contact impulses (the physical force chain, for the
     *              fracture model and visualisation) and does NOT move any floe.
     *  - Dynamics: the (C) heuristic — frozen floes are treated as fixed rest anchors (inverse mass 0,
     *              velocity 0) so the mobile floes' solved velocities are consistent with what actually
     *              moves (kills the first-order interpenetration / dt collapse). COMMITS velocities; does
     *              NOT record impulses (those come from the Forces pass). */
    enum class Mode { Forces = 0, Dynamics = 1 };

    /*! OPTIMJAM warm-start: seed each contact's impulses from the previous step's solution for the same
     *  floe-pair (nearest contact-point position). A held jam barely changes between steps, so projected
     *  Gauss-Seidel resumes near-converged and reaches a low residual in a few sweeps instead of the
     *  thousands a cold start needs. begin_step() MUST be called once per time step (it rotates the
     *  previous/next caches, separately per mode); set_warm_start(false) falls back to the cold start. */
    void set_warm_start(bool v) { m_warm_start = v; }
    bool warm_start() const { return m_warm_start; }
    //! Total number of per-floe kinetic-energy clamps applied by the Dynamics pass (run stats).
    long energy_clamps() const { return m_energy_clamps; }
    void begin_step() {
        for (int mo = 0; mo < 2; ++mo) { m_warm_prev[mo].swap(m_warm_next[mo]); m_warm_next[mo].clear(); }
    }

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
    bool solve(TGraph const& graph, Mode mode, int& out_iters, T& out_residual, T& out_vmax_free)
    {
        using floe::lcp::builder::GraphLCP;

        const bool dynamics = (mode == Mode::Dynamics);
        auto& warm_prev = m_warm_prev[(int)mode];
        auto& warm_next = m_warm_next[(int)mode];

        GraphLCP<T, TGraph> glcp(graph); // assembles M, invM, J, D, mu (no big LCP matrix built)
        const int n = glcp.nb_floes;
        const int m = glcp.nb_contacts;
        out_iters = 0;
        out_residual = 0;
        out_vmax_free = 0;
        if (m == 0) return true;

        // Which floes are FIXED anchors in this solve. Obstacles always are. In the Dynamics pass, frozen
        // floes also are — the (C) heuristic: the freeze is decided BEFORE the solve (see LCP_manager), so
        // the mobile floes come out consistent with what actually moves (no first-order interpenetration
        // from a mover approaching a pinned neighbour the solve wrongly assumed was mobile). In the Forces
        // pass, frozen floes keep their TRUE mass and free velocity, so the physical force chain (load
        // pressing down to the boundary) is computed — only real obstacles are fixed there.
        std::vector<char> is_fixed(n, 0);
        for (int i = 0; i < n; ++i)
            is_fixed[i] = (graph[i].floe->is_obstacle() || (dynamics && graph[i].floe->state().is_jammed())) ? 1 : 0;

        // Free velocity v = current floe velocities (already include the wind/ocean drag of the
        // previous move; solve_contacts runs before the next move). This is the "restart" velocity.
        // Fixed floes are clamped to rest so the contact solve sees them as immovable anchors.
        ublas::vector<T> v(3 * n);
        for (int i = 0; i < n; ++i) {
            if (is_fixed[i]) { v(3 * i) = v(3 * i + 1) = v(3 * i + 2) = 0; continue; }
            auto const st = graph[i].floe->state();
            v(3 * i)     = st.speed.x;
            v(3 * i + 1) = st.speed.y;
            v(3 * i + 2) = st.rot;
        }
        ublas::vector<T> v_free_vec = v; // pre-resolution velocities, kept for the kinetic-energy guard below

        // Per-contact local data: the (up to) 6 dof rows (3 per floe), the normal and tangential
        // jacobian values, the inverse-mass on each row, and the diagonal Delassus blocks Wnn / Wtt.
        struct Row { int idx; T jn; T dtg; T im; };
        std::vector<std::array<Row, 6>> rows(m);
        std::vector<T> Wnn(m, 0), Wtt(m, 0), muv(m, 0);
        // Per-contact identity for warm-start matching across steps: the floe pair and the contact point.
        std::vector<WarmKey> ckey(m);
        std::vector<T> cx(m, 0), cy(m, 0);
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
                        // Fixed (obstacle or frozen) floes have zero inverse mass: impulses don't move them.
                        const T im  = is_fixed[idx / 3] ? T(0) : glcp.invM(idx, idx);
                        rows[a][k] = Row{ idx, jn, dtg, im };
                        Wnn[a] += im * jn * jn;
                        Wtt[a] += im * dtg * dtg;
                    }
                    muv[a] = glcp.mu(a, a);
                    auto const& cp = graph[edge][c];
                    ckey[a] = make_key(cp.floe1, cp.floe2);
                    auto const cc = cp.frame.center();
                    cx[a] = cc.x; cy[a] = cc.y;
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

        // OPTIMJAM warm-start: seed the contact impulses from the previous step's committed solution
        // (matched per floe-pair, then by nearest contact-point position), and bring v into consistency
        // with the seeded impulses so the first sweep's relative velocities J*v are correct (PGS invariant
        // v = v_free + W p). Computed AFTER vmax_free/tol_feasible above, which must use the *free* velocity.
        if (m_warm_start && !warm_prev.empty())
        {
            for (int a = 0; a < m; ++a)
            {
                auto it = warm_prev.find(ckey[a]);
                if (it == warm_prev.end()) continue;
                const WarmContact* best = nullptr; T best_d2 = std::numeric_limits<T>::max();
                for (auto const& w : it->second) {
                    const T d2 = (w.x - cx[a]) * (w.x - cx[a]) + (w.y - cy[a]) * (w.y - cy[a]);
                    if (d2 < best_d2) { best_d2 = d2; best = &w; }
                }
                if (best) { pn(a) = best->pn; pt(a) = best->pt; }
            }
            for (int a = 0; a < m; ++a)
                if (pn(a) != T(0) || pt(a) != T(0))
                    for (int k = 0; k < 6; ++k) {
                        auto const& r = rows[a][k];
                        v(r.idx) += r.im * (r.jn * pn(a) + r.dtg * pt(a));
                    }
        }

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

        // OPTIMJAM energy guard (Dynamics pass only): an inelastic (e=0) contact solve must NOT inject
        // kinetic energy. A non-converged sweep or a divergent impulse (a floe squeezed between near-rigid
        // anchors) can make a single floe "explode" and shake/collapse the pack. Cap any floe whose post-
        // solve kinetic energy exceeds the WHOLE component's free (pre-solve) kinetic energy — no single
        // floe can carry more energy than the system started with. This is surgical: gross outliers are
        // scaled back, normal floes (and legitimate momentum transfer below that bound) are untouched. It
        // mirrors the kinetic-energy non-increase condition the Lemke path enforces (LCP_solver, calcEc>1).
        // KE is proportional to v^2 / invM (mass = 1/invM); the common 1/2 factor cancels in the comparison.
        if (dynamics)
        {
            T ke_free_tot = 0;
            for (int i = 0; i < n; ++i) if (!is_fixed[i])
                for (int d = 0; d < 3; ++d) {
                    const T iM = glcp.invM(3*i+d, 3*i+d);
                    if (iM > 0) ke_free_tot += v_free_vec(3*i+d) * v_free_vec(3*i+d) / iM;
                }
            if (ke_free_tot > 0)
                for (int i = 0; i < n; ++i) if (!is_fixed[i])
                {
                    T ke_i = 0;
                    for (int d = 0; d < 3; ++d) {
                        const T iM = glcp.invM(3*i+d, 3*i+d);
                        if (iM > 0) ke_i += v_best(3*i+d) * v_best(3*i+d) / iM;
                    }
                    if (ke_i > ke_free_tot) {
                        const T s = std::sqrt(ke_free_tot / ke_i);
                        v_best(3*i) *= s; v_best(3*i+1) *= s; v_best(3*i+2) *= s;
                        ++m_energy_clamps; // stats: reported at end of run
                    }
                }
        }

        // DYNAMICS pass only: commit the best velocities (this is what the move uses). Obstacles keep their
        // state; frozen floes are set to rest (held immobile in the anchored structure — this also stops the
        // wind/ocean drag from accumulating a phantom velocity in a floe that never moves, which would
        // otherwise fire it out on release). Only mobile floes receive their solved velocity. The FORCES pass
        // does NOT touch velocities, so the Forces pass (run first) reads the intact free velocities.
        if (dynamics)
            for (int i = 0; i < n; ++i)
            {
                if (graph[i].floe->is_obstacle()) continue;
                auto& st = graph[i].floe->state();
                if (st.is_jammed()) { st.speed = { 0, 0 }; st.rot = 0; continue; }
                st.speed = { v_best(3 * i), v_best(3 * i + 1) };
                st.rot   = v_best(3 * i + 2);
            }

        // FORCES pass only: feed the physical contact impulses back into the graph (the force chain that the
        // fracture model and the visualisation read). The Dynamics pass must NOT do this — its (C) infinite
        // masses produce non-physical, sometimes divergent impulses (a floe squeezed between rigid anchors).
        if (!dynamics)
        {
            ublas::vector<T> normal(m), tangential(2 * m);
            for (int a = 0; a < m; ++a) {
                normal(a)          = pn_best(a);
                tangential(2 * a)     = pt_best(a);
                tangential(2 * a + 1) = 0;
            }
            glcp.apply_impulses(normal, tangential);
        }

        // OPTIMJAM warm-start: record this pass's impulses for the next step's seed (per-mode cache).
        if (m_warm_start)
            for (int a = 0; a < m; ++a)
                warm_next[ckey[a]].push_back(WarmContact{ cx[a], cy[a], pn_best(a), pt_best(a) });

        return true; // committed (best-effort)
    }

private:
    int m_max_iter; //!< maximum number of sweeps
    T   m_tol;      //!< convergence tolerance (velocity-scaled impulse change of a sweep)

    // OPTIMJAM warm-start state. A contact is identified across steps by its (unordered) floe pointer pair;
    // a pair may carry several contact points, disambiguated by position at lookup. Caches are rotated by
    // begin_step() once per time step: solve() reads warm_prev and appends to warm_next, so the several
    // components solved within one step do not clobber each other (their keys are disjoint). One pair of
    // caches PER MODE (Forces / Dynamics), since the two passes have different (finite- vs infinite-mass)
    // impulse solutions — indexed by (int)Mode.
    bool m_warm_start{true};                                       //!< master switch for warm-starting
    long m_energy_clamps{0};                                       //!< total per-floe energy clamps (stats)
    struct WarmContact { T x, y, pn, pt; };                        //!< stored contact point + impulses
    using WarmKey = std::pair<const void*, const void*>;           //!< unordered floe-pointer pair
    std::map<WarmKey, std::vector<WarmContact>> m_warm_prev[2], m_warm_next[2];
    static WarmKey make_key(const void* a, const void* b) { return (a < b) ? WarmKey{a, b} : WarmKey{b, a}; }
};

}}} // namespace floe::lcp::solver

#endif // FLOE_LCP_SOLVER_GS_SOLVER_HPP
