/*!
 * \file ope/LCP_manager.hpp
 * \brief LCP manager
 * \author Quentin Jouet
 */

#ifndef OPE_LCP_MANAGER_HPP
#define OPE_LCP_MANAGER_HPP

#include "floe/ope/LCP_solver.hpp"
// #include <boost/numeric/ublas/vector_proxy.hpp> //remove?
// #include "floe/lcp/builder/graph_to_lcp.hpp" //remove?
#include "floe/ope/time_scale_manager.hpp"
// #include <boost/graph/graph_utility.hpp> //remove?

#include <iostream> // debug

namespace floe { namespace ope
{

/*! LCPManager
 *
 * Operator for LCP processing
 *
 */


class LCPManager
{

public:
    using solver_type = LCPSolver;
    using real = double; // TODO do better

    inline solver_type& get_solver() { return m_solver; }

    template<typename TContactGraph>
    void solve_contacts(TContactGraph& contact_graph);

private:

    solver_type m_solver;

    // template<typename TGraph, typename TGraphLCP, typename TLCP>
    // void update_floes_state(TGraph& graph, const TGraphLCP& graph_lcp, TLCP& lcp);
    template<typename TGraph>
    void update_floes_state(TGraph& graph, const boost::numeric::ublas::vector<real> Sol);
};


template<typename TContactGraph>
void LCPManager::solve_contacts(TContactGraph& contact_graph)
{
    auto const subgraphs = collision_subgraphs( contact_graph );
    int LCP_count = 0, nb_success = 0;
    bool r = 1;
    for ( auto& subgraph : subgraphs )
    {
        auto asubgraphs = active_subgraphs( subgraph );
        int loop_cnt = 0;
        while (asubgraphs.size() != 0 && r && loop_cnt < 100)
        {
            LCP_count += asubgraphs.size();
            for ( auto const& graph : asubgraphs )
            {
            //     floe::lcp::builder::GraphLCP<real, decltype(graph)> graph_lcp( graph );
            //     auto lcp = graph_lcp.getLCP();
            //     r = m_solver.solve( lcp );
                bool success;
                auto Sol = m_solver.solve( graph, success );
                if (success) nb_success++;
                update_floes_state(graph, Sol);
            }
            asubgraphs = active_subgraphs( subgraph );
            loop_cnt++;
        }
    }
    if (LCP_count)
        std::cout << " #LCP solve: "<< nb_success << " / " <<LCP_count <<std::endl;
}


// template<typename TGraph, typename TGraphLCP, typename TLCP>
// void LCPManager::update_floes_state(TGraph& graph, const TGraphLCP& graph_lcp, TLCP& lcp){

//     using namespace boost::numeric::ublas;

//     std::size_t m = graph_lcp.J.size2(); // num contacts (todo : cleaner way)
//     auto Sol = prod(
//         graph_lcp.invM,
//         prod(graph_lcp.J, subrange(lcp.z, 0, m)) + prod(graph_lcp.D, subrange(lcp.z, m, 3*m))
//     );

//     // std::size_t n = graph_lcp.J.size1() / 3; // num floes (todo : cleaner way)
//     for ( auto const v : boost::make_iterator_range( vertices(graph) ) )
//     {
//         graph[v]->state().speed += {Sol(3*v), Sol(3*v + 1)};
//         graph[v]->state().rot += Sol(3*v + 2);
//     }
// }

template<typename TGraph>
void LCPManager::update_floes_state(TGraph& graph, const boost::numeric::ublas::vector<real> Sol){

    for ( auto const v : boost::make_iterator_range( vertices(graph) ) )
    {
        graph[v]->state().speed = {Sol(3*v), Sol(3*v + 1)};
        graph[v]->state().rot = Sol(3*v + 2);
    }
}



}} // namespace floe::ope


#endif // OPE_LCP_MANAGER_HPP