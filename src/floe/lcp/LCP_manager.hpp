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
    inline void set_gs_params(int min_contacts, int max_iter, real_type rel_speed_max, real_type freeze_disp) {
        if (min_contacts > 0)  m_gs_min_contacts = min_contacts;
        m_gs_solver.set_max_iter(max_iter);
        if (rel_speed_max > 0) m_gs_rel_speed_max = rel_speed_max;
        if (freeze_disp > 0)   m_gs_freeze_disp = freeze_disp;
        if (m_optim_jam)
            std::cout << "OPTIMJAM params: min_contacts=" << m_gs_min_contacts
                      << " gs_max_iter=" << m_gs_solver.max_iter()
                      << " rel_speed_max=" << m_gs_rel_speed_max
                      << " freeze_disp=" << m_gs_freeze_disp << std::endl;
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
    real_type m_gs_freeze_disp{1e-3};             //!< freeze a floe whose predicted step displacement |v|*dt is below this fraction of its diameter
    solver::GaussSeidelSolver<real_type> m_gs_solver; //!< the alternative contact solver

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
            int gs_iters = 0;
            real_type gs_resid = 0, gs_vmax = 0;
            if ( m_gs_solver.solve( subgraph, gs_iters, gs_resid, gs_vmax ) )
            {
                long frozen = 0;
                // Freeze the MOVE only of floes that would not actually advance this step: those whose
                // predicted displacement |v|*dt is below a small fraction of their own diameter. This is
                // the physically meaningful "cannot move" test (a cantilever floe with real tangential
                // velocity advances and is NOT frozen, so it rolls/falls into place), and it self-regulates
                // the Zeno collapse: when dt shrinks because some floe is geometrically stuck, |v|*dt drops
                // below the threshold for everyone and the pack freezes for that step (breaking the Zeno),
                // then dt recovers and the genuinely-moving floes advance again next step.
                (void)gs_vmax;
                if ( m_gs_freeze )
                {
                    for ( auto v : boost::make_iterator_range(vertices(subgraph)) )
                    {
                        auto* floe = subgraph[v].floe;
                        if ( floe->is_obstacle() ) continue;
                        const auto& sp = floe->state().speed;
                        const real_type step_disp = std::sqrt(sp.x*sp.x + sp.y*sp.y) * dt;
                        if ( step_disp < m_gs_freeze_disp * floe->min_diameter() )
                        {
                            floe->state().set_jammed(true);
                            ++frozen;
                        }
                    }
                }
                std::cout << "OPTIMJAM GS | floes=" << num_vertices(subgraph)
                          << " contacts=" << num_contacts(subgraph)
                          << " sweeps=" << gs_iters << " residual=" << gs_resid
                          << " frozen=" << frozen << "/" << num_vertices(subgraph)
                          << " (committed)" << std::endl;
                continue;
            }
            std::cout << "OPTIMJAM GS | floes=" << num_vertices(subgraph)
                      << " contacts=" << num_contacts(subgraph)
                      << " sweeps=" << gs_iters << " residual=" << gs_resid
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
