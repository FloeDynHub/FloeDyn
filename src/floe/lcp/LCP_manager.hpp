/*!
 * \file lcp/LCP_manager.hpp
 * \brief LCP manager
 * \author Quentin Jouet
 */

#ifndef OPE_LCP_MANAGER_HPP
#define OPE_LCP_MANAGER_HPP

#include "floe/domain/time_scale_manager.hpp"
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/blas.hpp>
#include <iostream> // debug
#include <array>
#include <algorithm>
#include <chrono> // test perf

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
            // double m_JW = (m_nb_lcp_failed_stats[4] == 0)? 0 : m_nb_lcp_failed_stats[5]/m_nb_lcp_failed_stats[4];
            // real_type m_Ccoul = (m_nb_lcp_failed_stats[6] == 0)? 0 : m_nb_lcp_failed_stats[7]/m_nb_lcp_failed_stats[6];
            std::cout << "#TOTAL LCP solve: " << m_nb_lcp_success << "/" << m_nb_lcp << "(" << success_ratio() << "%) \n"
                << "LCP_failed compression phase: " << m_nb_lcp_failed_stats[0] << ", LCP_failed decompression phase:" << m_nb_lcp_failed_stats[1] << "\n"
                << "Sol to keep EC+<=EC-: " << m_nb_lcp_failed_stats[2] << "\n"
                << "#LCP with ECc>1: " << m_nb_lcp_failed_stats[3] << "\n"
                << "#LCP with J^TW<0: " << m_nb_lcp_failed_stats[4] << "\n";
            if (m_nb_lcp_failed_stats[4]!=0){std::cout << "mean of J^TW<0: " << m_nb_lcp_failed_stats[5]/m_nb_lcp_failed_stats[4] << "\n";}    
                std::cout << "#LCP with Cond Coul<0: " << m_nb_lcp_failed_stats[6] << "\n";
            if (m_nb_lcp_failed_stats[6]!=0){std::cout << "mean of Ccoul<0: " << m_nb_lcp_failed_stats[7]/m_nb_lcp_failed_stats[6] << "\n";}
            if (m_nb_lcp_failed_stats[0]!=0){std::cout << "mean of Err: " << m_nb_lcp_failed_stats[8]/m_nb_lcp_failed_stats[0] << "\n";}

                // if (m_nb_lcp_success!=0){
                //     std::cout << "Mean of Condition Number of Delassus matrix for LCP solved: " << m_nb_lcp_failed_stats[3]/(m_nb_lcp_success+m_nb_lcp_failed_stats[1]-m_nb_lcp_failed_stats[6]) << "\n"
                //     << "critical case:" << m_nb_lcp_failed_stats[6] << "\n";
                // }
                // if (m_nb_lcp_failed_stats[0]!=0){
                //     std::cout << "Mean of Condition Number of Delassus matrix for LCP failed during the compression phase: " << m_nb_lcp_failed_stats[4]/(m_nb_lcp_failed_stats[0]-m_nb_lcp_failed_stats[7]) << "\n"
                //     << "critical case:" << m_nb_lcp_failed_stats[7] << "\n";
                // }
                // if (m_nb_lcp_failed_stats[1]!=0){
                //     std::cout << "Mean of Condition Number of Delassus matrix for LCP failed during the decompression phase: " << m_nb_lcp_failed_stats[5]/(m_nb_lcp_failed_stats[1]-m_nb_lcp_failed_stats[8]) << "\n"
                //     << "critical case:" << m_nb_lcp_failed_stats[8] << "\n";
                // }  
            // printf ("#TOTAL LCP solve : %ld / %ld (%f %%) \n", m_nb_lcp_success, m_nb_lcp, success_ratio());
                // LCP_failed compression phase: %d, decompression phase: %d \n 
                // Sol to keep EC+<=EC-: %d \n", m_nb_lcp_success, m_nb_lcp, success_ratio(), m_nb_lcp_failed_comp, m_nb_lcp_failed_decomp, m_nb_lcp_E_sup );
    }

    //! LCP solver accessor
    inline solver_type& get_solver() { return m_solver; }

    //! Solve collision represented by a contact graph
    template<typename TContactGraph>
    int solve_contacts(TContactGraph& contact_graph);
    //! Get solving success ratio in percent
    double success_ratio(){ return (m_nb_lcp == 0)? 100 : 100 * (double)m_nb_lcp_success/m_nb_lcp; }

private:

    solver_type m_solver; //!< LCP Solver
    long m_nb_lcp; //!< Total number of LCP managed
    long m_nb_lcp_success; //!< Total number of LCP solving success
    //Matt
    real_type m_nb_lcp_failed_stats[9]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; // LCP failed statistics: [nb LCP failed during compression phase,
    // nb LCP failed during decompression phase, nb LCP using linear solution to keep EC+ <= EC-]
    // long m_nb_lcp_failed_decomp{0}; // total number of LCP failed during the decompression phase
    // long m_nb_lcp_E_sup{0}; // total number of LCP using linear solution to keep EC+ <= EC-
    //EndMatt

    double chrono_active_subgraph{0.0}; // test perf
    double max_chrono_active_subgraph{0.0}; // test perf

    //! Update floes state with LCP solution
    template<typename TContactGraph>
    void update_floes_state(TContactGraph& graph, const std::array<value_vector, 2> Sol);
};


template<typename T>
template<typename TContactGraph>
int LCPManager<T>::solve_contacts(TContactGraph& contact_graph)
{
    auto const subgraphs = collision_subgraphs( contact_graph );
    int LCP_count = 0, nb_success = 0;
    double nb_lcp_failed_stats[9];
    for (int i=0;i<9;++i){
        nb_lcp_failed_stats[i] = 0;
    }
    
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
        

        // Active subhraph LCP strategy
        auto asubgraphs = active_subgraphs( subgraph );
        std::size_t loop_cnt = 0;
        int loop_nb_success = -1;
        while (asubgraphs.size() != 0
               && loop_cnt < std::min( 60 * num_contacts(subgraph), std::size_t{1000})
               && loop_nb_success !=0 )
        {
            int loop_nb_success = 0;
            LCP_count += asubgraphs.size();
            for ( auto const& graph : asubgraphs )
            {
                bool success;
                if (num_contacts(graph) > 50){
                    std::cout << " Q4," << std::flush;
                    for ( auto const& igraph : quad_cut( graph ) ){
                        auto Sol = m_solver.solve( igraph, success, nb_lcp_failed_stats );
                        mark_solved(igraph, success);
                        if (success) loop_nb_success++;
                        update_floes_state(igraph, Sol);
                    }
                } else {
                    auto Sol = m_solver.solve( graph, success, nb_lcp_failed_stats );
                    mark_solved(graph, success);
                    if (success) loop_nb_success++;
                    update_floes_state(graph, Sol);
                }
                mark_changed_parent(graph, subgraph);
            }
            // auto t_start2 = std::chrono::high_resolution_clock::now(); // test perf
            asubgraphs = active_subgraphs( subgraph );
            // auto t_end2 = std::chrono::high_resolution_clock::now(); // test perf
            // auto call_time = std::chrono::duration<double, std::milli>(t_end2-t_start2).count(); // test perf
            // chrono_active_subgraph += call_time; // test perf
            // max_chrono_active_subgraph = std::max(max_chrono_active_subgraph, call_time); // test perf
            nb_success += loop_nb_success;
            loop_cnt++;
        }
        if (asubgraphs.size() != 0)
        {
            LCP_count += asubgraphs.size();
            for ( auto const& graph : asubgraphs ) mark_solved(graph, false);
        }

        // {
        //     // Big LCP solving
        //     LCP_count += 1;
        //     bool success;
        //     auto Sol = m_solver.solve( subgraph, success );
        //     mark_solved(subgraph, success);
        //     std::cout << "BIG LCP ";
        //     if (success) {nb_success++; std::cout << "SUCCESS\n" ; }
        //     update_floes_state(subgraph, Sol);
        // }

    //     nb_active_subgraph_loop += loop_cnt; // test perf
    //     std::cout << "NB_LOOP : " << loop_cnt << " / " << 60 * num_contacts(subgraph) << "\n"; // test perf
    }
    // auto t_end = std::chrono::high_resolution_clock::now(); // test perf
    // auto call_time = std::chrono::duration<double, std::milli>(t_end-t_start).count(); // test perf
    // std::cout << "#slv() : " << m_solver.nb_solver_run // test perf
    // << ", T : " << (double)call_time // test perf
    // << ", avg_T_Slv : " << (double)m_solver.chrono_solver/m_solver.nb_solver_run // test perf
    // << ", max_T_Slv : " << m_solver.max_chrono_solver // test perf
    // << ", #ASG_loop : " << nb_active_subgraph_loop // test perf
    // << ", avg_T_ASG : " << chrono_active_subgraph/nb_active_subgraph_loop // test perf
    // << ", max_T_ASG : " << max_chrono_active_subgraph // test perf
    // << " ( #contacts : " << num_contacts(contact_graph) << " )" // test perf
    // << "\n"; // test perf

    // version omp
    // #pragma omp parallel for
    // for ( std::size_t i = 0; i < subgraphs.size(); ++i )
    // {
    //     auto& subgraph = subgraphs[i];
    //     auto asubgraphs = active_subgraphs( subgraph );
    //     int loop_cnt = 0;
    //     while (asubgraphs.size() != 0 && loop_cnt < 60 * num_contacts(subgraph) )
    //     {
    //         LCP_count += asubgraphs.size();
    //         #pragma omp parallel for
    //         for ( std::size_t j = 0; j < asubgraphs.size(); ++j )
    //         {
    //             bool success;
    //             auto& graph = asubgraphs[j];
    //             auto Sol = m_solver.solve( graph, success );
    //             mark_solved(graph, success);
    //             // #pragma omp critical
    //             if (success) 
    //             {
    //                 // #pragma omp atomic
    //                 // nb_success++;
    //                 update_floes_state(graph, Sol);
    //             }
    //         }
    //         asubgraphs = active_subgraphs( subgraph );
    //         loop_cnt++;
    //     }
    // }
    m_nb_lcp += LCP_count;
    m_nb_lcp_success += nb_success;
    for (int i=0;i<9;++i){
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
void LCPManager<T>::update_floes_state(TContactGraph& graph, const std::array<value_vector, 2> Sol){

    for ( auto const v : boost::make_iterator_range( vertices(graph) ) )
    {
        graph[v].floe->state().speed = {Sol[0](3*v), Sol[0](3*v + 1)}; // fv_test
        graph[v].floe->state().rot = Sol[0](3*v + 2); // fv_test
        graph[v].floe->add_impulse(Sol[1](v)); // fv_test
    }
}



}} // namespace floe::lcp


#endif // OPE_LCP_MANAGER_HPP