/*!
 * \file ope/LCP_manager.hpp
 * \brief LCP manager
 * \author Quentin Jouet
 */

#ifndef OPE_LCP_MANAGER_HPP
#define OPE_LCP_MANAGER_HPP

#include "floe/ope/LCP_solver.hpp"
#include "floe/ope/time_scale_manager.hpp"

#include <iostream> // debug

#ifdef _OPENMP
#include <omp.h>
#endif

namespace floe { namespace ope
{

/*! LCPManager
 *
 * Operator for LCP processing
 *
 */


template<typename T>
class LCPManager
{

public:
    using value_type = T;
    using solver_type = LCPSolver<value_type>;

    ~LCPManager(){ printf ("#TOTAL LCP solve : %d / %d (%f %%) \n", m_nb_lcp_success, m_nb_lcp, success_ratio() ); }

    inline solver_type& get_solver() { return m_solver; }

    template<typename TContactGraph>
    void solve_contacts(TContactGraph& contact_graph);
    double success_ratio(){ return (m_nb_lcp == 0)? 100 : 100 * (double)m_nb_lcp_success/m_nb_lcp; }

private:

    solver_type m_solver;
    int m_nb_lcp;
    int m_nb_lcp_success;

    template<typename TContactGraph>
    void update_floes_state(TContactGraph& graph, const boost::numeric::ublas::vector<value_type> Sol);
};


template<typename T>
template<typename TContactGraph>
void LCPManager<T>::solve_contacts(TContactGraph& contact_graph)
{
    auto const subgraphs = collision_subgraphs( contact_graph );
    int LCP_count = 0, nb_success = 0;
    for ( auto& subgraph : subgraphs )
    {
        auto asubgraphs = active_subgraphs( subgraph );
        int loop_cnt = 0;
        while (asubgraphs.size() != 0 && loop_cnt < 60 * num_contacts(subgraph) )
        {
            LCP_count += asubgraphs.size();
            for ( auto const& graph : asubgraphs )
            {
                bool success;
                auto Sol = m_solver.solve( graph, success );
                mark_solved(graph, success);
                if (success) 
                {
                    nb_success++;
                    update_floes_state(graph, Sol);
                }
            }
            asubgraphs = active_subgraphs( subgraph );
            loop_cnt++;
        }
    }
    /* version omp
    #pragma omp parallel for
    for ( std::size_t i = 0; i < subgraphs.size(); ++i )
    {
        auto& subgraph = subgraphs[i];
        auto asubgraphs = active_subgraphs( subgraph );
        int loop_cnt = 0;
        while (asubgraphs.size() != 0 && loop_cnt < 60 * num_contacts(subgraph) )
        {
            LCP_count += asubgraphs.size();
            // pragma omp parallel for
            for ( std::size_t j = 0; j < asubgraphs.size(); ++j )
            {
                bool success;
                auto& graph = asubgraphs[j];
                auto Sol = m_solver.solve( graph, success );
                mark_solved(graph, success);
                #pragma omp critical
                if (success) 
                {
                    nb_success++;
                    update_floes_state(graph, Sol);
                }
            }
            asubgraphs = active_subgraphs( subgraph );
            loop_cnt++;
        }
    }*/
    m_nb_lcp += LCP_count;
    m_nb_lcp_success += nb_success;
    if (LCP_count)
        std::cout << " #LCP solve: "<< nb_success << " / " << LCP_count << std::endl;
}


template<typename T>
template<typename TContactGraph>
void LCPManager<T>::update_floes_state(TContactGraph& graph, const boost::numeric::ublas::vector<value_type> Sol){

    for ( auto const v : boost::make_iterator_range( vertices(graph) ) )
    {
        graph[v]->state().speed = {Sol(3*v), Sol(3*v + 1)};
        graph[v]->state().rot = Sol(3*v + 2);
    }
}



}} // namespace floe::ope


#endif // OPE_LCP_MANAGER_HPP