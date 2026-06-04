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

#include <boost/graph/adjacency_list.hpp>

// saving matrix when lcp solver failed for further analysing
#include "H5Cpp.h"
#include "floe/utils/random.hpp"                    // for saving lcp statistic matrix


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
    }

    //! LCP solver accessor
    inline solver_type& get_solver() { return m_solver; }

    //! OPTIMJAM: enable/disable the jamming optimization (off => behaviour identical to baseline)
    inline void set_optim_jam(bool v) {
        m_optim_jam = v;
        std::cout << "OPTIMJAM optimization is " << (v ? "ENABLED" : "disabled") << std::endl;
    }
    inline bool optim_jam() const { return m_optim_jam; }

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

    bool   m_optim_jam{false}; //!< OPTIMJAM: master switch for the jamming optimization
    long   m_jam_frozen_total{0}; //!< OPTIMJAM: number of floes frozen on the previous step (for transition prints)

    //! OPTIMJAM: detect jammed clusters and freeze them (see implementation for the full rationale).
    template<typename TSubgraphs>
    void detect_and_freeze_jam(TSubgraphs& subgraphs, real_type dt);

    //! Update floes state with LCP solution
    template<typename TContactGraph>
    void update_floes_state(TContactGraph& graph, const value_vector Sol, real_type time);
    //! Update floes impulses with LCP solution
    template<typename TContactGraph>
    void update_floes_impulses(TContactGraph& graph, real_type time, real_type dt);

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


/*! OPTIMJAM — jamming detection & freezing.
 *
 *  Called once per step, BEFORE the LCP loop. It decides, floe by floe, which floes are "jammed"
 *  and freezes them (velocity = 0, jammed flag set). A jammed floe is then treated as a fixed wall
 *  (infinite mass) by the LCP builder (see graph_to_lcp.hpp), so:
 *    - a fully jammed cluster has no active contact left  => its LCP is skipped entirely;
 *    - in a hybrid cluster (floes still raining on a jammed pack), only the moving region and its
 *      interface with the pack are solved, the deep frozen core drops out.
 *
 *  Detection is per-floe (NOT per-subgraph-hash): a new floe landing on the pack does not reset the
 *  state of the floes already jammed underneath. A floe is frozen when, for JAM_CONFIRM_STEPS
 *  consecutive steps, it is slow (kinetic energy per unit mass below JAM_V2_THRESHOLD) AND anchored
 *  (its connected component contains an obstacle — a force chain can only hold against a fixed boundary).
 *
 *  "Restart before resolve" is automatic: every step PROBLEM::step_solve() calls unjam_all_floes()
 *  and the move step applies wind/ocean drag to every floe (jammed ones keep their position frozen but
 *  still accumulate the drag in their velocity). Hence whenever the LCP does run, the floes carry their
 *  correct free velocity (the wind/current push) and the solver produces the true load-bearing forces
 *  against the obstacles — never the spurious "a floe hits a group at rest" problem.
 *
 *  Un-jamming: every JAM_PROBE_PERIOD steps a confirmed-jammed floe is left UNfrozen for one step, so
 *  the LCP re-evaluates it with the current load. If it comes back moving (e.g. an arch that yielded
 *  under accumulated pressure), its counter resets and it is released; otherwise it re-freezes.
 *
 *  \warning v1 validity domain: meant for spatially & temporally uniform forcing, fracture OFF.
 *           Spontaneous (noise-triggered) collapse of a cluster under constant load+composition is
 *           suppressed by the freeze; collapse here is load-triggered (re-evaluated on the probe step).
 */
template<typename T>
template<typename TSubgraphs>
void LCPManager<T>::detect_and_freeze_jam(TSubgraphs& subgraphs, real_type /*dt*/)
{
    const real_type JAM_V2_THRESHOLD  = 1e-2; //!< kinetic_energy/mass threshold (~|v| < 0.14 m/s) for "slow"
    const int       JAM_CONFIRM_STEPS = 20;   //!< consecutive slow+anchored steps before freezing (cf. the "20x" memory)
    const int       JAM_PROBE_PERIOD  = 100;  //!< re-evaluate a frozen floe with a real LCP solve every N steps

    long frozen_now = 0;

    for (auto& subgraph : subgraphs)
    {
        // A cluster can only jam if it is anchored to a fixed boundary (obstacle).
        bool anchored = false;
        for (auto v : boost::make_iterator_range(vertices(subgraph)))
            if (subgraph[v].floe->is_obstacle()) { anchored = true; break; }

        for (auto v : boost::make_iterator_range(vertices(subgraph)))
        {
            auto* floe = subgraph[v].floe;
            if (floe->is_obstacle()) continue;
            auto& st = floe->state();

            const real_type v2 = floe->kinetic_energy() / floe->mass();
            if (anchored && v2 < JAM_V2_THRESHOLD)
                st.inc_jam_count();
            else
                st.reset_jam_count();

            const bool confirmed = anchored && st.get_jam_count() >= JAM_CONFIRM_STEPS;
            const bool probe     = confirmed && (st.get_jam_count() % JAM_PROBE_PERIOD == 0);
            if (confirmed && !probe)
            {
                st.speed = {0, 0};
                st.rot   = 0;
                st.set_jammed(true);
                ++frozen_now;
            }
        }
    }

    if (frozen_now != m_jam_frozen_total)
    {
        std::cout << "OPTIMJAM | frozen floes: " << m_jam_frozen_total << " -> " << frozen_now
                  << (frozen_now > m_jam_frozen_total ? "  (jam growing)" : "  (jam releasing)")
                  << std::endl;
    }
    m_jam_frozen_total = frozen_now;
}


template<typename T>
template<typename TContactGraph>
int LCPManager<T>::solve_contacts(TContactGraph& contact_graph, real_type time, real_type dt)
{

    auto const subgraphs = collision_subgraphs( contact_graph );
    int LCP_count=0, nb_success=0;
    int nb_lcp_failed_stats[3]={0,0,0};

    const std::size_t limit_sup_loop_cnt    = 800;//5000; // from Quentin: 1000
    const std::size_t limit_sup_nb_contact  =  800;//500; // from Quentin:   50

    // OPTIMJAM: decide & freeze jammed floes before solving. Freezing makes them infinite-mass walls
    // for the LCP, so fully-frozen clusters generate no active subgraph (their LCP is skipped) and
    // hybrid clusters only solve their moving region. No effect when the optimization is disabled.
    if (m_optim_jam)
        detect_and_freeze_jam( subgraphs, dt );

    // variables for contact informations:
    #ifdef LCPSTATS
        static bool end_recording = false;
    #endif

    for ( auto& subgraph : subgraphs )
    {
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
    update_floes_impulses(contact_graph, time, dt);
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
void LCPManager<T>::update_floes_impulses(TContactGraph& graph, real_type time, real_type dt){
    for ( auto const v : boost::make_iterator_range( vertices(graph) ) )
    {
        auto* floe = graph[v].floe;
        const real_type imp = graph[v].impulse(); // fv_test

        if (m_optim_jam && floe->state().is_jammed() && imp == 0)
        {
            // Deep-frozen floe that was excluded from the LCP this step: re-inject the cached
            // impulse rate so its cumulative impulse ("red") keeps growing at the correct rate.
            floe->add_impulse(floe->get_jam_impulse_rate() * dt);
        }
        else
        {
            floe->add_impulse(imp);
            // Remember the per-second impulse while the floe is actually solved (and not frozen),
            // so we have a fresh value to re-inject once it freezes.
            if (m_optim_jam && !floe->state().is_jammed() && dt > 0)
                floe->set_jam_impulse_rate(imp / dt);
        }
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
