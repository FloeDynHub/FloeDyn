/*!
 * \file floe/lcp/jam_manager.hpp
 * \brief OPTIMJAM — dense-jam handling: routing, temporal freeze & probes, two-pass Gauss-Seidel
 *        solve, guard-rails and run statistics.
 * \author Quentin Jouet
 *
 * Extracted from LCPManager (pure move, behaviour identical) so that solve_contacts keeps reading
 * as the historical Lemke resolution, with a single OPTIMJAM entry point: try_solve_component().
 *
 * Full reference (motivation, pipeline, parameters, error budget): doc/optimjam.md. In short:
 *  - large, ANCHORED, QUASI-STATIC contact components (dense jams) are routed to a projected
 *    Gauss-Seidel solve of the whole component (degeneracy-tolerant, see GS_solver.hpp) instead of
 *    the Lemke + active-subgraph cascade (problem A: LCP count explosion);
 *  - floes detected as geometrically blocked — real velocity but ~no NET displacement over several
 *    steps — are FROZEN (their position integration is skipped), which stops the drag-driven creep
 *    that collapses the adaptive time step (problem B: Zeno); a staggered periodic PROBE releases
 *    them to re-test their mobility (that is how piles/arches collapse);
 *  - the freeze is decided BEFORE the solve and frozen floes are treated as fixed anchors in the
 *    DYNAMICS pass (velocities consistent with what actually moves), while a second FORCES pass at
 *    true masses on the free velocities records the physical force chain (fracture model + viz);
 *  - guard-rails in the FloeDyn spirit (never block): best-effort commits, Lemke fallback on solver
 *    decline, per-floe kinetic-energy clamp, anti-limit-cycle emergency freeze on repeated
 *    no-progress state recoveries.
 */

#ifndef OPE_JAM_MANAGER_HPP
#define OPE_JAM_MANAGER_HPP

#include "../product/config/config_base.hpp" // types
#include "floe/lcp/solver/GS_solver.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/range/iterator_range.hpp>

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <limits>

namespace floe { namespace lcp
{

using namespace types;

class JamManager
{

public:

    /*! Master switch (CLI --optim_jam). When disabled, try_solve_component() always declines and
     *  the LCPManager behaves exactly like the baseline. */
    inline void set_enabled(bool v) {
        m_enabled = v;
        std::cout << "OPTIMJAM (Gauss-Seidel path) is " << (v ? "ENABLED" : "disabled") << std::endl;
    }
    inline bool enabled() const { return m_enabled; }

    //! Routing/solver/freeze tuning (CLI --jam_params).
    //! \param min_contacts   minimum number of contacts for a component to be routed to Gauss-Seidel
    //! \param max_iter       maximum number of Gauss-Seidel sweeps
    //! \param rel_speed_max  route only components whose max contact approach speed (m/s) is below this
    //! \param eps            net-progress threshold, fraction of the floe diameter (temporal freeze)
    //! \param stuck_N        consecutive no-progress steps before freezing a floe
    //! \param probe_K        probe period: every K steps release a frozen floe to retest its mobility
    inline void set_params(int min_contacts, int max_iter, real_type rel_speed_max,
                           real_type eps, int stuck_N, int probe_K) {
        if (min_contacts > 0)  m_gs_min_contacts = min_contacts;
        m_gs_solver.set_max_iter(max_iter);
        if (rel_speed_max > 0) m_gs_rel_speed_max = rel_speed_max;
        if (eps > 0)           m_gs_eps = eps;
        if (stuck_N > 0)       m_gs_stuck_N = stuck_N;
        if (probe_K > 0)       m_gs_probe_K = probe_K;
        if (m_enabled)
            std::cout << "OPTIMJAM params: min_contacts=" << m_gs_min_contacts
                      << " gs_max_iter=" << m_gs_solver.max_iter()
                      << " rel_speed_max=" << m_gs_rel_speed_max
                      << " eps=" << m_gs_eps
                      << " stuck_N=" << m_gs_stuck_N
                      << " probe_K=" << m_gs_probe_K << std::endl;
    }

    /*! Freeze the *move* of floes the temporal criterion detects as blocked (CLI --jam_freeze).
     *
     *  The GS solver gives the correct held contact forces but does not stop the floes from slowly
     *  creeping closer under the wind/ocean drag applied during the move (drag is added after the
     *  contact solve). Over time that creep closes the gaps, which (a) collapses the adaptive time
     *  step (Zeno) and (b) produces exactly-overlapping contacts with a degenerate (NaN) normal.
     *  Skipping the position integration of blocked floes (via the is_jammed flag, re-decided every
     *  step) prevents the creep, hence both pathologies. Disable to compare with GS alone. */
    inline void set_freeze(bool v) { m_gs_freeze = v; }

    /*! Enable/disable warm-starting the Gauss-Seidel solver across time steps (default on).
     *  Seeds each contact's impulses from the previous step's solution so a held jam — whose
     *  configuration barely changes step to step — resumes near-converged and needs far fewer
     *  sweeps. See GS_solver.hpp. */
    inline void set_warm_start(bool v) { m_gs_solver.set_warm_start(v); }

    /*! Guard-rail: called by the problem each time "dt too small" forces a state recovery
     *  (safe_move_floe_group / detect_proximity). If recovery repeats WITHOUT simulated-time
     *  progress, the run is trapped in a deterministic INTER/RECOVER limit cycle (observed in
     *  practice: with --bustle 0 everything is deterministic, the recovery restores the exact same
     *  state, and the only evolving state — the per-floe stuck counters — drives a probe pattern
     *  that is periodic mod probe_K). Response, in the same spirit as the active-subgraph loop
     *  guard (exit unresolved so life goes on): freeze the floes already showing a blockage for the
     *  next few steps, so the jam stops re-trying the moves that fail, dt recovers and time
     *  advances; then resume normal probing/freezing. Imperfect but unblocking. */
    inline void notify_recover(real_type time) {
        if (!m_enabled) return;
        if (time <= m_last_recover_time) {  // recovered again without advancing past the previous recovery
            if (++m_recover_repeat >= 2 && m_emergency_steps_left == 0) {
                m_emergency_steps_left = 20;
                ++m_nb_emergencies;
                std::cout << "OPTIMJAM EMERGENCY: recover loop detected (no time progress) -> "
                          << "freezing the stuck floes (counter>=1) of routed components for "
                          << m_emergency_steps_left << " steps" << std::endl;
            }
        } else {
            m_recover_repeat = 0;
        }
        m_last_recover_time = time;
    }

    /*! Once per time step (start of solve_contacts): rotate the Gauss-Seidel warm-start caches
     *  (last step's impulses become this step's seed) and tick down the anti-limit-cycle emergency
     *  window (see notify_recover). */
    inline void begin_step() {
        if (!m_enabled) return;
        m_gs_solver.begin_step();
        if (m_emergency_steps_left > 0) --m_emergency_steps_left;
    }

    /*! The single OPTIMJAM entry point, called by LCPManager::solve_contacts for each connected
     *  contact component. Returns true if the component was handled here (freeze decided, dynamics
     *  pass committed, force chain recorded) — the caller skips it. Returns false (component left
     *  untouched and unfrozen) when OPTIMJAM is disabled, the component is not routed (small, free
     *  or non-quasi-static: impulsive collisions stay on Lemke), or the solver declined — the
     *  caller then runs the historical Lemke resolution. */
    template<typename TSubgraph>
    bool try_solve_component(TSubgraph const& subgraph);

    //! End-of-run summary (printed next to the historical #TOTAL LCP block).
    ~JamManager() {
        // These numbers are the measurable counterpart of the tuning knobs (jam_params) for the
        // fidelity/speed trade-off.
        if (m_gs_routed)
            std::cout   << "#OPTIMJAM GS: " << m_gs_committed << "/" << m_gs_routed << " components committed, "
                        << m_gs_fallback << " Lemke fallbacks\n"
                        << "  dynamics pass: avg sweeps " << (m_gs_committed ? m_gs_dyn_sweeps / m_gs_committed : 0)
                        << ", saturated (best-effort commit) " << m_gs_dyn_saturated << "/" << m_gs_committed << "\n"
                        << "  forces pass:   avg sweeps " << (m_gs_committed ? m_gs_frc_sweeps / m_gs_committed : 0)
                        << ", saturated (chain not fully converged) " << m_gs_frc_saturated << "/" << m_gs_committed << "\n"
                        << "  energy clamps (floes): " << m_gs_solver.energy_clamps()
                        << " | emergency en-bloc freezes: " << m_nb_emergencies << "\n";
    }

private:

    //! Should this contact component be routed to the Gauss-Seidel solver ?
    //! Route only large, anchored (touching an obstacle) and quasi-static components: that is exactly
    //! the dense-jam regime where Lemke + active-subgraph explodes and fails to converge.
    template<typename TSubgraph>
    bool route_to_gs(TSubgraph const& subgraph) const;

    bool      m_enabled{false};                   //!< master switch for the GS path
    int       m_gs_min_contacts{50};              //!< route a component to GS above this many contacts
    real_type m_gs_rel_speed_max{0.5};            //!< route only components with max contact approach speed below this (m/s)
    bool      m_gs_freeze{true};                  //!< freeze the move of blocked floes (prevents creep/Zeno/NaN)
    // Temporal jam-detection knobs (replace the failed per-instant |v|*dt criterion): freeze a floe whose
    // net displacement stays below m_gs_eps*diameter for m_gs_stuck_N consecutive steps; every m_gs_probe_K
    // steps skip the freeze (staggered per floe) to let a released floe prove it can move again.
    real_type m_gs_eps{3e-4};                     //!< net-progress threshold, as a fraction of the floe diameter
    int       m_gs_stuck_N{10};                   //!< consecutive no-progress steps before freezing
    int       m_gs_probe_K{10};                   //!< probe period: every K steps, release to retest mobility
    solver::GaussSeidelSolver<real_type> m_gs_solver; //!< the alternative contact solver

    // Anti-limit-cycle guard-rail state (see notify_recover)
    int       m_emergency_steps_left{0};          //!< >0: freeze the stuck floes of routed components for that many steps
    int       m_recover_repeat{0};                //!< consecutive recoveries without simulated-time progress
    real_type m_last_recover_time{-std::numeric_limits<real_type>::max()}; //!< sim time of the last recovery
    long      m_nb_emergencies{0};                //!< stats: how many times the guard-rail fired

    // End-of-run stats (printed by the destructor, like the historical #TOTAL LCP block)
    long m_gs_routed{0};         //!< components routed to the GS path
    long m_gs_committed{0};      //!< components the GS dynamics pass committed
    long m_gs_fallback{0};       //!< components declined by GS -> Lemke fallback
    long m_gs_dyn_sweeps{0};     //!< cumulated dynamics-pass sweeps (avg = /m_gs_committed)
    long m_gs_frc_sweeps{0};     //!< cumulated forces-pass sweeps
    long m_gs_dyn_saturated{0};  //!< dynamics passes that hit the sweep cap (best-effort commits)
    long m_gs_frc_saturated{0};  //!< forces passes that hit the sweep cap (recorded chain not fully converged)
};


template<typename TSubgraph>
bool JamManager::route_to_gs(TSubgraph const& subgraph) const
{
    if ((int)num_contacts(subgraph) < m_gs_min_contacts) return false;

    // Must be anchored to a fixed boundary (obstacle): a force chain can only hold against one.
    bool anchored = false;
    for (auto v : boost::make_iterator_range(vertices(subgraph)))
        if (subgraph[v].floe->is_obstacle()) { anchored = true; break; }
    if (!anchored) return false;

    // Route only if the contacts are quasi-static: max contact *approach* speed below threshold. This
    // gates on relative contact velocity (held / slowly grinding) rather than absolute floe velocity,
    // so an anchored pack that is barely rearranging routes to GS even if its kinetic energy is high.
    real_type max_approach = 0;
    for (auto e : boost::make_iterator_range(edges(subgraph)))
        for (std::size_t i = 0; i < subgraph[e].size(); ++i)
        {
            const real_type rs = subgraph[e][i].relative_speed(); // < 0 means approaching
            if (-rs > max_approach) max_approach = -rs;
        }
    return max_approach < m_gs_rel_speed_max;
}


template<typename TSubgraph>
bool JamManager::try_solve_component(TSubgraph const& subgraph)
{
    if ( !m_enabled || !route_to_gs( subgraph ) ) return false;

    const bool gs_emergency = ( m_emergency_steps_left > 0 );

    // Decide the per-floe freeze BEFORE the solve, so the Gauss-Seidel solve is consistent with
    // what actually moves (it treats frozen floes as fixed anchors — see GS_solver). This removes
    // the dominant, first-order interpenetration the old solve-then-freeze order caused: a mobile
    // floe approaching a neighbour the solve assumed was mobile but which we then pinned.
    //
    // The freeze criterion is TEMPORAL, not per-instant. A jammed floe has a real (drift) velocity —
    // it *wants* to move — but is held by its neighbours, so no instantaneous velocity/|v|*dt test
    // can tell "moving" from "wants-to-move-but-blocked". The reliable signal is the ACTUAL NET
    // DISPLACEMENT over several steps: a floe whose net progress stays below eps*diameter for stuck_N
    // consecutive steps is blocked -> freeze it. A frozen floe makes 0 progress, so it would stay
    // stuck forever; hence every probe_K steps we skip the freeze (staggered per floe by vertex index,
    // to avoid a synchronous mass release that would re-jam) to let a floe whose constraint was
    // released show it can move again — and the consistent solve then sorts it out: a still-blocked
    // probe gets v~0 (no INTER), a probe with free space actually moves (the pile can collapse).
    ++m_gs_routed;
    long frozen = 0;
    long dbg_reset = 0; int dbg_maxcnt = 0; real_type dbg_maxratio = 0, dbg_diam = 0; // DIAG
    if ( m_gs_freeze )
    {
        for ( auto v : boost::make_iterator_range(vertices(subgraph)) )
        {
            auto* floe = subgraph[v].floe;
            if ( floe->is_obstacle() ) continue;
            // Anti-limit-cycle emergency (see notify_recover): freeze the floes already showing
            // a blockage (stuck counter >= 1), probes suppressed for them, counters paused. The
            // storm culprits (wedged floes, e.g. against a wall) have high counters -> covered.
            // Floes flowing freely (counter 0 — e.g. a lattice of floes falling onto the jam,
            // which belongs to the same routed component) are NOT pinned mid-flight: en-bloc
            // freezing them caused visible artifacts (suspended blocks, shear gaps opening in
            // the falling lattice), since a frozen floe is also committed at rest (v=0).
            if ( gs_emergency && floe->jam_tracked() && floe->jam_stuck_counter() >= 1 )
                { floe->state().set_jammed(true); ++frozen; continue; }
            dbg_diam = floe->static_floe().max_diameter(); // characteristic floe size (min_diameter() returns ~0)
            const auto cur = floe->state().real_position();
            if ( !floe->jam_tracked() ) // first observation: capture the reference, don't freeze
            {
                floe->set_jam_ref_pos(cur);
                floe->set_jam_tracked(true);
                floe->set_jam_stuck_counter(0);
                continue;
            }
            const auto ref = floe->jam_ref_pos();
            const real_type dx = cur.x - ref.x, dy = cur.y - ref.y;
            const real_type net_disp = std::sqrt(dx*dx + dy*dy);
            const real_type thresh = m_gs_eps * dbg_diam; // = m_gs_eps * max_diameter (already computed above)
            if ( thresh > 0 && net_disp / thresh > dbg_maxratio ) dbg_maxratio = net_disp / thresh; // DIAG
            if ( net_disp > thresh ) // genuine net progress -> moving
            {
                floe->set_jam_ref_pos(cur);
                floe->set_jam_stuck_counter(0);
                ++dbg_reset; // DIAG
                continue;
            }
            const int cnt = floe->jam_stuck_counter() + 1; // no net progress this step
            floe->set_jam_stuck_counter(cnt);
            if ( cnt > dbg_maxcnt ) dbg_maxcnt = cnt; // DIAG
            const bool probe = ( (cnt + (int)v) % m_gs_probe_K == 0 );
            if ( cnt >= m_gs_stuck_N && !probe )
            {
                floe->state().set_jammed(true);
                ++frozen;
            }
        }
    }

    // (beta) TWO passes per component (see GS_solver::Mode):
    //  - DYNAMICS first: the (C) heuristic (frozen floes = fixed rest anchors) gives the velocities
    //    we actually move with -> stable dt.
    //  - FORCES second: a TRUE-mass solve on the floes' FREE (pre-resolution) velocities computes the
    //    physical force chain we record (fracture + visualisation), without moving anyone. The force a
    //    held floe transmits comes from its free, drag-loaded velocity, so the force pass must read
    //    those — but the dynamics pass overwrites the velocities, so we snapshot/restore around it.
    // The force pass runs only when the dynamics pass committed; otherwise the component falls back to
    // Lemke, which records its own impulses (so no double counting), and the freezes are undone.
    using gs_mode = typename decltype(m_gs_solver)::Mode;

    std::vector<std::array<real_type,3>> v_free; // free (pre-resolution) velocities of the component
    v_free.reserve(num_vertices(subgraph));
    for ( auto v : boost::make_iterator_range(vertices(subgraph)) ) {
        auto const& st = subgraph[v].floe->state();
        v_free.push_back({ st.speed.x, st.speed.y, st.rot });
    }

    int d_it = 0; real_type d_resid = 0, d_vmax = 0;
    if ( m_gs_solver.solve( subgraph, gs_mode::Dynamics, d_it, d_resid, d_vmax ) )
    {
        (void)d_vmax;
        ++m_gs_committed;
        m_gs_dyn_sweeps += d_it;
        if ( d_it >= m_gs_solver.max_iter() ) ++m_gs_dyn_saturated;
        // Capture the committed dynamics velocities, restore the free velocities for the force pass,
        // run the force pass (records impulses, moves nobody), then restore the dynamics velocities.
        std::vector<std::array<real_type,3>> v_dyn;
        v_dyn.reserve(num_vertices(subgraph));
        std::size_t k = 0;
        for ( auto v : boost::make_iterator_range(vertices(subgraph)) ) {
            auto& st = subgraph[v].floe->state();
            v_dyn.push_back({ st.speed.x, st.speed.y, st.rot });
            st.speed = { v_free[k][0], v_free[k][1] }; st.rot = v_free[k][2]; ++k;
        }
        int f_it = 0; real_type f_resid = 0, f_vmax = 0;
        m_gs_solver.solve( subgraph, gs_mode::Forces, f_it, f_resid, f_vmax ); // records the force chain
        m_gs_frc_sweeps += f_it;
        if ( f_it >= m_gs_solver.max_iter() ) ++m_gs_frc_saturated;
        k = 0;
        for ( auto v : boost::make_iterator_range(vertices(subgraph)) ) {
            auto& st = subgraph[v].floe->state();
            st.speed = { v_dyn[k][0], v_dyn[k][1] }; st.rot = v_dyn[k][2]; ++k;
        }
        std::cout << "OPTIMJAM GS | floes=" << num_vertices(subgraph)
                  << " contacts=" << num_contacts(subgraph)
                  << " dyn_sweeps=" << d_it << " dyn_resid=" << d_resid
                  << " frc_sweeps=" << f_it << " frc_resid=" << f_resid
                  << " frozen=" << frozen << "/" << num_vertices(subgraph)
                  << " | reset=" << dbg_reset << " maxcnt=" << dbg_maxcnt
                  << " maxratio=" << dbg_maxratio << " diam=" << dbg_diam // DIAG
                  << " (committed)" << std::endl;
        return true;
    }
    // Dynamics pass declined: undo the freezes we set for this component so the Lemke fallback runs its
    // baseline (it solves and moves all the floes, and records its own impulses). Counters stay as set.
    ++m_gs_fallback;
    for ( auto v : boost::make_iterator_range(vertices(subgraph)) )
        subgraph[v].floe->state().set_jammed(false);
    std::cout << "OPTIMJAM GS | floes=" << num_vertices(subgraph)
              << " contacts=" << num_contacts(subgraph)
              << " dyn_sweeps=" << d_it << " dyn_resid=" << d_resid
              << " (NOT converged -> Lemke fallback)" << std::endl;
    return false; // fall through to the Lemke active-subgraph resolution
}

}} // namespace floe::lcp


#endif // OPE_JAM_MANAGER_HPP
