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

// #include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/graph_utility.hpp>
// #include <boost/graph/filtered_graph.hpp>
// #include <boost/graph/connected_components.hpp>
// #include <boost/graph/copy.hpp>

// #include "floe/collision/contact_point.hpp"
// #include "floe/collision/floe_contact.hpp"
// #include "floe/collision/floe_vertex.hpp"


#ifdef _OPENMP
#include <omp.h>
#endif

namespace floe { namespace lcp
{

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
    using real_type = typename solver_type::real_type;
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

    //! Solve collision represented by a contact graph
    template<typename TContactGraph>
    int solve_contacts(TContactGraph& contact_graph, real_type time);
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
template<typename TContactGraph>
int LCPManager<T>::solve_contacts(TContactGraph& contact_graph, typename T::real_type time)
{

    auto const subgraphs = collision_subgraphs( contact_graph );
    int LCP_count=0, nb_success=0;
    int nb_lcp_failed_stats[3]={0,0,0}; 

    const std::size_t limit_sup_loop_cnt    = 800;//5000; // from Quentin: 1000
    const std::size_t limit_sup_nb_contact  =  80;//500; // from Quentin:   50

    // variables for contact informations:
    #ifdef LCPSTATS
        static bool end_recording = false;
    #endif

    // m_solver.nb_solver_run = 0; // test perf
    // m_solver.chrono_solver = 0; // test perf
    // m_solver.max_chrono_solver = 0; // test perf
    // chrono_active_subgraph = 0; // test perf
    // max_chrono_active_subgraph = 0; // test perf
    // int nb_active_subgraph_loop = 0;// test
    // auto t_start = std::chrono::high_resolution_clock::now(); // test perf
    for ( auto& subgraph : subgraphs )
    {
        //  // Big LCP solving
        // LCP_count += 1;
        // bool success;
        // auto Sol = m_solver.solve( subgraph, success );
        // mark_solved(subgraph, success);
        // if (success) nb_success++;
        // update_floes_state(subgraph, Sol, subgraph);

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
            // auto t_start2 = std::chrono::high_resolution_clock::now(); // test perf
            asubgraphs = active_subgraphs( subgraph ); // computes the new relative velocitoies from velocities of modified floes 
            
            // auto t_end2 = std::chrono::high_resolution_clock::now(); // test perf
            // auto call_time = std::chrono::duration<double, std::milli>(t_end2-t_start2).count(); // test perf
            // chrono_active_subgraph += call_time; // test perf
            // max_chrono_active_subgraph = std::max(max_chrono_active_subgraph, call_time); // test perf
            nb_success += loop_nb_success;
            ++loop_cnt;

            if (loop_nb_success==0) {contact_loop_stats[1] = 0;}
        }
        if (asubgraphs.size() != 0)
        {
            all_solved = false;
            if (loop_nb_success!=0) {contact_loop_stats[1] = 2;}
            std::cout << "End of the while loop without resolution of all contacts!! nb contact: "<< num_contacts(subgraph) << "\n";
            // LCP_count += asubgraphs.size();

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