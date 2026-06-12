/*!
 * \file lcp/LCP_manager.hpp
 * \brief LCP manager
 * \author Quentin Jouet
 */

#ifndef OPE_LCP_MANAGER_HPP
#define OPE_LCP_MANAGER_HPP

#include "floe/domain/time_scale_manager.hpp"
#include "../product/config/config_base.hpp" // types

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/blas.hpp>
#include <iostream> // debug
#include <array>
#include <algorithm>
#include <chrono> // test perf
#include <vector>
#include <limits>

#include <boost/graph/adjacency_list.hpp>

// saving matrix when lcp solver failed for further analysing
#include "H5Cpp.h"
#include "floe/utils/random.hpp"                    // for saving lcp statistic matrix

#include "floe/lcp/solver/GS_solver.hpp"            // OPTIMJAM: alternative Gauss-Seidel contact solver


#ifdef _OPENMP
#include <omp.h>
#endif

namespace floe { namespace lcp
{

using namespace types;

/*! LCPManager
 *
 * Operator for Collision processing
 *
 */


template<typename TSolver>
class LCPManager
{

public:
    using solver_type = TSolver;
    using value_vector = boost::numeric::ublas::vector<real_type>;

    //! Constructor
    LCPManager(real_type epsilon) : m_solver{epsilon}, m_nb_lcp{0}, m_nb_lcp_success{0} {}

    //! Destructor
    ~LCPManager(){
        if (m_nb_lcp)
            std::cout   << "#TOTAL LCP solve: " << m_nb_lcp_success << "/" << m_nb_lcp << "(" << success_ratio() << "%) \n"
                        << "LCP_failed compression phase: " << m_nb_lcp_failed_stats[0] <<
                        ", LCP_failed decompression phase: " << m_nb_lcp_failed_stats[1] << "\n" <<
                        ", LCP_solved with solution maintaining the kinetic energy: " << m_nb_lcp_failed_stats[2] << "\n";
        // OPTIMJAM end-of-run summary (mirrors the #TOTAL LCP block above): how often the GS path ran, how
        // healthy its two passes were, and how often the guard-rails had to intervene. These numbers are the
        // measurable counterpart of the tuning knobs (jam_params) for the fidelity/speed trade-off.
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

    //! LCP solver accessor
    inline solver_type& get_solver() { return m_solver; }

    /*! OPTIMJAM: enable/disable the alternative Gauss-Seidel path.
     *
     *  When enabled, large *anchored* and *quasi-static* contact components (dense jams) are routed to
     *  a single projected Gauss-Seidel solve (see GS_solver.hpp) instead of the Lemke + active-subgraph
     *  resolution. Everything else (impulsive collisions) keeps using the historical Lemke solver.
     *  When disabled, behaviour is identical to baseline.
     */
    inline void set_optim_jam(bool v) {
        m_optim_jam = v;
        std::cout << "OPTIMJAM (Gauss-Seidel path) is " << (v ? "ENABLED" : "disabled") << std::endl;
    }
    inline bool optim_jam() const { return m_optim_jam; }

    //! OPTIMJAM: routing/solver tuning.
    //! \param min_contacts   minimum number of contacts for a component to be routed to Gauss-Seidel
    //! \param max_iter       maximum number of Gauss-Seidel sweeps
    //! \param rel_speed_max  route only components whose max contact approach speed (m/s) is below this
    //! \param eps            net-progress threshold, fraction of the floe diameter (temporal freeze)
    //! \param stuck_N        consecutive no-progress steps before freezing a floe
    //! \param probe_K        probe period: every K steps release a frozen floe to retest its mobility
    inline void set_gs_params(int min_contacts, int max_iter, real_type rel_speed_max,
                              real_type eps, int stuck_N, int probe_K) {
        if (min_contacts > 0)  m_gs_min_contacts = min_contacts;
        m_gs_solver.set_max_iter(max_iter);
        if (rel_speed_max > 0) m_gs_rel_speed_max = rel_speed_max;
        if (eps > 0)           m_gs_eps = eps;
        if (stuck_N > 0)       m_gs_stuck_N = stuck_N;
        if (probe_K > 0)       m_gs_probe_K = probe_K;
        if (m_optim_jam)
            std::cout << "OPTIMJAM params: min_contacts=" << m_gs_min_contacts
                      << " gs_max_iter=" << m_gs_solver.max_iter()
                      << " rel_speed_max=" << m_gs_rel_speed_max
                      << " eps=" << m_gs_eps
                      << " stuck_N=" << m_gs_stuck_N
                      << " probe_K=" << m_gs_probe_K << std::endl;
    }

    /*! OPTIMJAM: freeze the *move* of floes of a component the GS solver resolved as a held equilibrium.
     *
     *  The GS solver gives the correct held contact forces but does not stop the floes from slowly
     *  creeping closer under the wind/ocean drag applied during the move (drag is added after the contact
     *  solve). Over time that creep closes the gaps, which (a) collapses the adaptive time step (Zeno) and
     *  (b) produces exactly-overlapping contacts with a degenerate (NaN) normal. Skipping the position
     *  integration of a GS-confirmed held component (via the existing is_jammed flag, re-decided every
     *  step) prevents the creep, hence both pathologies. Disable to compare with GS alone. */
    inline void set_gs_freeze(bool v) { m_gs_freeze = v; }

    /*! OPTIMJAM: enable/disable warm-starting the Gauss-Seidel solver across time steps (default on).
     *  Seeds each contact's impulses from the previous step's solution so a held jam — whose configuration
     *  barely changes step to step — resumes near-converged and needs far fewer sweeps. See GS_solver.hpp. */
    inline void set_gs_warm_start(bool v) { m_gs_solver.set_warm_start(v); }

    /*! OPTIMJAM guard-rail: called by the problem each time "dt too small" forces a state recovery
     *  (safe_move_floe_group / detect_proximity). If recovery repeats WITHOUT simulated-time progress, the
     *  run is trapped in a deterministic INTER/RECOVER limit cycle (observed in practice: with --bustle 0
     *  everything is deterministic, the recovery restores the exact same state, and the only evolving state
     *  — the per-floe stuck counters — drives a probe pattern that is periodic mod probe_K). Response, in
     *  the same spirit as the active-subgraph loop guard (exit unresolved so life goes on): freeze the
     *  routed components EN BLOC (suspend probes) for the next few steps, so nothing in the jam moves, dt
     *  recovers and time advances; then resume normal probing/freezing. Imperfect but unblocking. */
    inline void notify_recover(real_type time) {
        if (!m_optim_jam) return;
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

    //! Solve collision represented by a contact graph
    template<typename TContactGraph>
    int solve_contacts(TContactGraph& contact_graph, real_type time, real_type dt);
    //! Get solving success ratio in percent
    double success_ratio(){ return (m_nb_lcp == 0)? 100 : 100 * (double)m_nb_lcp_success/m_nb_lcp; }

private:

    solver_type m_solver; //!< LCP Solver
    long        m_nb_lcp; //!< Total number of LCP managed
    long        m_nb_lcp_success; //!< Total number of LCP solving success
    long        m_nb_lcp_failed_stats[3]={0,0,0}; // LCP failed statistics: [nb LCP failed during compression phase,
    // nb LCP failed during decompression phase, nb LCP solved maintaining the kinetic energy in decompression phase].
    // Other LCP statistics are found in the matlab routine (see folder: io/outputs).

    double chrono_active_subgraph{0.0}; // test perf
    double max_chrono_active_subgraph{0.0}; // test perf

    // OPTIMJAM — Gauss-Seidel alternative solver and routing
    bool      m_optim_jam{false};                 //!< master switch for the GS path
    int       m_gs_min_contacts{50};              //!< route a component to GS above this many contacts
    real_type m_gs_rel_speed_max{0.5};            //!< route only components with max contact approach speed below this (m/s)
    bool      m_gs_freeze{true};                  //!< freeze the move of GS-confirmed held floes (prevents creep/Zeno/NaN)
    // Temporal jam-detection knobs (replace the failed per-instant |v|*dt criterion): freeze a floe whose
    // net displacement stays below m_gs_eps*diameter for m_gs_stuck_N consecutive steps; every m_gs_probe_K
    // steps skip the freeze (staggered per floe) to let a released floe prove it can move again.
    real_type m_gs_eps{3e-4};                     //!< net-progress threshold, as a fraction of the floe diameter
    int       m_gs_stuck_N{10};                   //!< consecutive no-progress steps before freezing
    int       m_gs_probe_K{10};                   //!< probe period: every K steps, release to retest mobility
    solver::GaussSeidelSolver<real_type> m_gs_solver; //!< the alternative contact solver

    // OPTIMJAM anti-limit-cycle guard-rail state (see notify_recover)
    int       m_emergency_steps_left{0};          //!< >0: freeze routed components en bloc for that many steps
    int       m_recover_repeat{0};                //!< consecutive recoveries without simulated-time progress
    real_type m_last_recover_time{-std::numeric_limits<real_type>::max()}; //!< sim time of the last recovery
    long      m_nb_emergencies{0};                //!< stats: how many times the guard-rail fired

    // OPTIMJAM end-of-run stats (printed by the destructor, like the historical #TOTAL LCP block)
    long m_gs_routed{0};         //!< components routed to the GS path
    long m_gs_committed{0};      //!< components the GS dynamics pass committed
    long m_gs_fallback{0};       //!< components declined by GS -> Lemke fallback
    long m_gs_dyn_sweeps{0};     //!< cumulated dynamics-pass sweeps (avg = /m_gs_committed)
    long m_gs_frc_sweeps{0};     //!< cumulated forces-pass sweeps
    long m_gs_dyn_saturated{0};  //!< dynamics passes that hit the sweep cap (best-effort commits)
    long m_gs_frc_saturated{0};  //!< forces passes that hit the sweep cap (recorded chain not fully converged)

    //! OPTIMJAM: should this contact component be routed to the Gauss-Seidel solver ?
    //! Route only large, anchored (touching an obstacle) and quasi-static components: that is exactly
    //! the dense-jam regime where Lemke + active-subgraph explodes and fails to converge.
    template<typename TSubgraph>
    bool route_to_gs(TSubgraph const& subgraph) const;

    //! Update floes state with LCP solution
    template<typename TContactGraph>
    void update_floes_state(TContactGraph& graph, const value_vector Sol, real_type time);
    //! Update floes impulses with LCP solution
    template<typename TContactGraph>
    void update_floes_impulses(TContactGraph& graph, real_type time);

    /*! \fn bool saving_contact_graph_in_hdf5(int lCP_count, std::size_t loop_count, std::size_t size_a_sub_graph, bool all_solved )
        \brief Saves information on the contact graph in the same file as LCP statistics.

        \param lcp_count the total number of treated LCP at the end of the "while loop"
        \param loop_count the total number of loop performed at the end of the "while loop"
        \param size_asubgraph the number of LCP generated from the first active sub-graphic
        \param all_soved a boolean to indicate whether all treated LCP during the "while loop" are solved

        \details the data are saved in the <a>nx6</a> matrix form, with <a>n</a> the number of sub-graph dealt with during the simulation.
        The output is a boolean to ensure data do not exceed a defined maximum storage.
        Each line of the matrix contains: the numerous of the last unsolved LCP stored in the h5 file (see LCPsolver::saving_LCP_in_hdf5),
        the numerous of the last solved LCP stored in the h5 file, lcp_count, loop_count, size_asubgraph and all_solved.

        \warning Do not be confused with contact graph, sub-graph and active sub-graph (see LCPManager::solve_contacts).
    */
    bool saving_contact_graph_in_hdf5(int lcp_count, std::size_t loop_count, std::size_t size_asubgraph,
        bool all_solved, int contact_loop_stats[] );
};


template<typename T>
template<typename TSubgraph>
bool LCPManager<T>::route_to_gs(TSubgraph const& subgraph) const
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


template<typename T>
template<typename TContactGraph>
int LCPManager<T>::solve_contacts(TContactGraph& contact_graph, real_type time, real_type dt)
{

    auto const subgraphs = collision_subgraphs( contact_graph );
    int LCP_count=0, nb_success=0;
    int nb_lcp_failed_stats[3]={0,0,0};

    // OPTIMJAM: rotate the Gauss-Seidel warm-start cache once per step (last step's impulses become the
    // seed for this step's solves; the new step's solutions are accumulated for the next one), and tick
    // down the anti-limit-cycle emergency window (see notify_recover).
    if ( m_optim_jam ) {
        m_gs_solver.begin_step();
        if ( m_emergency_steps_left > 0 ) --m_emergency_steps_left;
    }
    const bool gs_emergency = ( m_emergency_steps_left > 0 );

    const std::size_t limit_sup_loop_cnt    = 800;//5000; // from Quentin: 1000
    const std::size_t limit_sup_nb_contact  =  800;//500; // from Quentin:   50

    // variables for contact informations:
    #ifdef LCPSTATS
        static bool end_recording = false;
    #endif

    for ( auto& subgraph : subgraphs )
    {
        // OPTIMJAM: route large, anchored, quasi-static components to the Gauss-Seidel solver, which
        // solves the whole component in one (degeneracy-tolerant) pass instead of thousands of Lemke
        // active-subgraph solves. The GS solver is a *safe fast-path*: it commits a solution only if it
        // reaches feasibility (no penetration velocity); otherwise it leaves the floes untouched and we
        // fall through to the historical Lemke resolution below. Everything else uses Lemke too.
        if ( m_optim_jam && route_to_gs( subgraph ) )
        {
            // OPTIMJAM (C): decide the per-floe freeze BEFORE the solve, so the Gauss-Seidel solve is
            // consistent with what actually moves (it treats frozen floes as fixed anchors — see GS_solver).
            // This removes the dominant, first-order interpenetration the old solve-then-freeze order caused:
            // a mobile floe approaching a neighbour the solve assumed was mobile but which we then pinned.
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
                continue;
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
            // fall through to the Lemke active-subgraph resolution below
        }

        // Active subgraph LCP strategy
        auto asubgraphs = active_subgraphs( subgraph );
        std::size_t loop_cnt    = 0;
        int loop_nb_success     = -1;
        bool active_quad_cut    = 0;

        // variables for contact informations:
        #ifdef LCPSTATS
            std::size_t size_a_sub_graph = asubgraphs.size();
        #endif
        bool all_solved = true;

        int contact_loop_stats[2]={0,0};    // number of contact points, indicator for be out of loop due to all success (1) or no success (0)
                                            // or no enough iteration (2)
        contact_loop_stats[0] = static_cast<int>(num_contacts(subgraph));
        contact_loop_stats[1] = 1;

        while (asubgraphs.size() != 0
               && loop_cnt < std::min( 100 * num_contacts(subgraph), limit_sup_loop_cnt) // 60 * num_contacts(subgraph)
               && loop_nb_success !=0 )
        {
            loop_nb_success = 0;    // if no succes after one total path of contact graph, no use (nothing change)
                                    // to browse again the while loop. Thus we add "loop_nb_succes !=0".

            LCP_count += asubgraphs.size();

            for ( auto const& graph : asubgraphs ) // loop over the total number of active contact group
            {
                bool success;
                if (num_contacts(graph) > limit_sup_nb_contact){
                    std::cout << "Q4, nb contact:" << num_contacts(graph) << " \n";
                    auto qdct = quad_cut( graph );
                    std::cout << "nb quad_cut: " << qdct.size() << " \n";
                    active_quad_cut = 1;
                    for ( auto const& igraph : quad_cut( graph ) ){
                        auto Sol = m_solver.solve( igraph, success, nb_lcp_failed_stats );
                        mark_solved(igraph, success);
                        if (success) {++loop_nb_success;}
                        update_floes_state(igraph, Sol, time);
                    }
                } else {
                    auto Sol = m_solver.solve( graph, success, nb_lcp_failed_stats );
                    mark_solved(graph, success);
                    if (success) {++loop_nb_success;}
                    update_floes_state(graph, Sol, time); // updates the velocity of floes
                }

                mark_changed_parent(graph, subgraph); // indicates which floes have been modified
            }
            asubgraphs = active_subgraphs( subgraph ); // computes the new relative velocitoies from velocities of modified floes
            nb_success += loop_nb_success;
            ++loop_cnt;

            if (loop_nb_success==0) {contact_loop_stats[1] = 0;}
        }
        if (asubgraphs.size() != 0)
        {
            all_solved = false;
            if (loop_nb_success!=0) {contact_loop_stats[1] = 2;}
            std::cout << "End of the while loop without resolution of all contacts!! nb contact: "<< num_contacts(subgraph) << "\n";

            for ( auto const& graph : asubgraphs ) mark_solved(graph, false);
        }

        // Mat
        // Saving data on LCP:
        /*
         * Recovery of contact data (LCP_count, etc). Save in h5 file:
         */
        #ifdef LCPSTATS
            if (!end_recording && size_a_sub_graph!=0) {
                end_recording = saving_contact_graph_in_hdf5( LCP_count, loop_cnt, size_a_sub_graph, all_solved, contact_loop_stats );
            }
        #endif
        // End saving data on LCP
        // EndMat
    }
    update_floes_impulses(contact_graph, time);
    m_nb_lcp += LCP_count;
    m_nb_lcp_success += nb_success;
    for (int i=0;i<3;++i){
        m_nb_lcp_failed_stats[i] += nb_lcp_failed_stats[i];
    }

    #ifndef MPIRUN
    if (LCP_count)
        std::cout << " #LCP solve: "<< nb_success << " / " << LCP_count << std::endl;
    #endif
    return nb_success;
}


template<typename T>
template<typename TContactGraph>
void LCPManager<T>::update_floes_state(TContactGraph& graph, const value_vector Sol, real_type time){

    for ( auto const v : boost::make_iterator_range( vertices(graph) ) )
    {
        graph[v].floe->state().speed = {Sol(3*v), Sol(3*v + 1)}; // fv_test
        graph[v].floe->state().rot = Sol(3*v + 2); // fv_test
    }
}

//! Update floes impulses from contact graph
template<typename T>
template<typename TContactGraph>
void LCPManager<T>::update_floes_impulses(TContactGraph& graph, real_type time){
    for ( auto const v : boost::make_iterator_range( vertices(graph) ) )
    {
        graph[v].floe->add_impulse(graph[v].impulse()); // fv_test
    }
    for ( auto const& edge : make_iterator_range( edges( graph ) ) )
    {
        for ( std::size_t i = 0; i < graph[edge].size(); ++i ) // iter over contacts
        {
            // Add contact point impulses to corresponding floes
            graph[source(edge, graph)].floe->add_contact_impulse(
                graph[edge][i].frame.center(), -graph[edge][i].impulse_abs_frame(),
                time);
            graph[target(edge, graph)].floe->add_contact_impulse(
                graph[edge][i].frame.center(), graph[edge][i].impulse_abs_frame(),
                time);
        }
    }
}


template<typename T>
bool LCPManager<T>::saving_contact_graph_in_hdf5(int LCP_count, std::size_t loop_count, std::size_t size_a_sub_graph,
 bool all_solved, int contact_loop_stats[] ){

    using namespace H5;                                         // proper way inside a function (to prevent extension to entire code)

    const H5std_string FILE_NAME("io/outputs/LCP_stats.h5");
    const H5std_string GROUP_NAME_I( "solved" ); // root group
    const H5std_string GROUP_NAME_II( "unsolved" ); // root group
    const H5std_string Last_Memb( "Last LCP" ); // to prevent similar LCP
    const H5std_string Contact_Graph_Info( "Contact Graph Info" );
    const hsize_t Max_storage_line = 50000;

    const int dim_stats = 8; // number of saved data on the contact loop.

    /*
     * Try block to detect exceptions raised by any of the calls inside it
     */
    try{
        /*
         * Turn off the auto-printing when failure occurs so that we can
         * handle the errors appropriately
         */
        Exception::dontPrint();
        /*
         * Create or Open a file.
         */
        H5File* file;
        try {
            file = new H5File( FILE_NAME, H5F_ACC_RDWR );
        } catch (...) {
            return false;
        }

        Group* M_solved = new Group(file->openGroup(GROUP_NAME_I));
        Group* M_unsolved = new Group(file->openGroup(GROUP_NAME_II));

        int last_lcp_uns[1];
        DataSet* dataset_LM_uns = new DataSet(M_unsolved->openDataSet( Last_Memb ));
        dataset_LM_uns->read( last_lcp_uns, PredType::NATIVE_INT );

        int last_lcp[1];
        DataSet* dataset_LM = new DataSet(M_solved->openDataSet( Last_Memb ));
        dataset_LM->read( last_lcp, PredType::NATIVE_INT );

        delete dataset_LM_uns;
        delete dataset_LM;

        DataSet* CGI;
        hsize_t dim_curr_cgi[2]={0,0};
        if (last_lcp_uns[0]==0 && last_lcp[0]==0) {
            delete M_unsolved;
            delete M_solved;
            delete file;
            return false;
        }
        else {
            try {
                CGI = new DataSet(file->openDataSet(Contact_Graph_Info));

                DataSpace fpsace_ind = CGI->getSpace();
                fpsace_ind.getSimpleExtentDims( dim_curr_cgi, NULL); // retrieves the current dimensions

                if (dim_curr_cgi[0] > Max_storage_line){
                    std::cout << "the maximum storage (" << Max_storage_line << ") for contact graph information is reached.\n";
                    delete M_unsolved;
                    delete M_solved;
                    delete CGI;
                    delete file;
                    return true;
                }

                hsize_t ext_size[2] = { dim_curr_cgi[0]+1, dim_curr_cgi[1] };
                CGI->extend( ext_size ); // extension with one new line

            }
            catch (...) {
                // Creation of dataset to store information on contact graph and the "while loop":
                hsize_t dim_cg[2] = {1, dim_stats};
                hsize_t maxdims_cg[2] = {H5S_UNLIMITED, dim_stats}; // unlimited dataspace
                DataSpace space_cg( 2, dim_cg, maxdims_cg );
                DSetCreatPropList prop_cg; // Modify dataset creation property to enable chunking
                hsize_t chunk_dims_cg[2] = {1, dim_stats}; // with extendible dataset we cannot use contiguous but chunked dataset
                prop_cg.setChunk(2, chunk_dims_cg);
                CGI = new DataSet(file->createDataSet( Contact_Graph_Info, PredType::NATIVE_INT, space_cg, prop_cg ));
            }
        }

        int contact_stat[dim_stats];
        contact_stat[0] = last_lcp_uns[0];
        contact_stat[1] = last_lcp[0];
        contact_stat[2] = LCP_count;
        contact_stat[3] = static_cast<int>(loop_count);
        contact_stat[4] = static_cast<int>(size_a_sub_graph);
        contact_stat[5] =  (all_solved == true)? 1:0;
        contact_stat[6] = contact_loop_stats[0]; // number of contact points
        contact_stat[7] = contact_loop_stats[1]; // indicator for be out of loop due to all success (1) or no success (0)
        // or no enough iteration (2)

        DataSpace fspace_cgi = CGI->getSpace();
        hsize_t dim_cgi[2] = {1,dim_stats};
        hsize_t offset_cgi[2] = {dim_curr_cgi[0], 0};
        fspace_cgi.selectHyperslab( H5S_SELECT_SET, dim_cgi, offset_cgi); // selection of the hyperslab
        DataSpace mspace_cgi( 2, dim_cgi );
        CGI->write(contact_stat, PredType::NATIVE_INT, mspace_cgi, fspace_cgi); // write in the hyperslab


        delete M_unsolved;
        delete M_solved;
        delete CGI;
        delete file;
    }  // end of try block
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
    error.printErrorStack();
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
    error.printErrorStack();
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
    error.printErrorStack();
    }
    return false;
}


}} // namespace floe::lcp


#endif // OPE_LCP_MANAGER_HPP
