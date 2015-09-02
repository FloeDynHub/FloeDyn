/*!
 * \file floe/collision/contact_graph.hpp
 * \brief Contact graph definition and manipulation.
 * \author Roland Denis
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/graph/doc/index.html">boost::graph</a>
 * \see <a href="http://stackoverflow.com/questions/26763193/return-a-list-of-connected-component-subgraphs-in-boost-graph">stackoverflow</a>
 */

#ifndef FLOE_COLLISION_CONTACT_GRAPH_HPP
#define FLOE_COLLISION_CONTACT_GRAPH_HPP

#include <cstddef>
#include <vector>
#include <functional>
#include <algorithm>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>

#include "floe/collision/contact_point.hpp"
#include "floe/collision/floe_contact.hpp"

#include <iostream> // DEBUG

namespace floe { namespace collision
{

using namespace boost;

/*! Contact Graph as an undirected adjacency_list with vertex and edge properties
 *
 * \tparam TContactPoint    Type of contact point.
 *
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/graph/doc/adjacency_list.html">boost/graph/adjacency_list.hpp</a>
 *
 * \remark Why not using std::list (or forward_list) for contact list ?
 */
template <
    typename TContactPoint
>
using ContactGraph = adjacency_list<
    vecS, // Out Edge List
    vecS, // Vertex List
    undirectedS, // (un)Directed
    typename TContactPoint::floe_type const*, // Vertex Property
    FloeContact<TContactPoint> // Edge Property
>;

/*! Return a filtered graph without obstacles
 *
 * \tparam  TGraph  Type of the graph.
 * \param   graph   The graph.
 * \return          A graph without obstacle vertices.
 */
template < typename TGraph >
inline
filtered_graph <
    TGraph, 
    boost::keep_all, 
    std::function< bool(typename TGraph::vertex_descriptor) > 
>
graph_without_obstacle( TGraph const& graph )
{
    return {
        graph,
        boost::keep_all(),
        [&graph] ( typename TGraph::vertex_descriptor v ) { return ! graph[v]->is_obstacle(); }
    };
}

/*! Return a filtered graph given a vertices mask ( mask[vertex] == value )
 *
 * \tparam  TGraph  Type of the graph.
 * \tparam  TMask   Type of the mask or value vector.
 * \tparam  T       Type of the value in the mask.
 * \param   graph   The graph.
 * \param   mask    The mask.
 * \param   value   The value that must match the mask.
 * \return          A graph with edge and vertices filtered.
 */
template < 
    typename TGraph,
    typename TMask,
    typename T
>
inline
filtered_graph <
    TGraph,
    std::function< bool( typename TGraph::edge_descriptor ) >,
    std::function< bool( typename TGraph::vertex_descriptor ) >
>
graph_from_mask( TGraph const& graph, TMask const& mask, T value = true  )
{
    return {
        graph,
        [&graph, mask, value] ( typename TGraph::edge_descriptor e ) { return mask[source(e, graph)] == value && mask[target(e, graph)] == value; },
        [mask, value] ( typename TGraph::vertex_descriptor v ) { return mask[v] == value; }
    };
}

/*! Returns a filtered graph given a list of vertices id
 *
 * \todo use boost::dynamic_bitset instead of vector<bool>
 * \tparam  TGraph  Type of the graph (auto-deduced).
 * \param   graph   The graph.
 * \param   ids     List of the vertices ids to filter.
 * \param   reverse The vertices are keeped if false, rejected if true.
 * \return          A graph with filtered edges and vertices.
 */
template <
    typename TGraph
>
inline
auto
graph_from_ids( TGraph const& graph, std::vector<std::size_t> const& ids, bool reverse = false )
    -> decltype( graph_from_mask( graph, std::vector<bool>(), true ) )
{
    std::vector<bool> mask( num_vertices(graph), false );
    for ( std::size_t i : ids )
        mask[i] = true;

    return graph_from_mask( graph, mask, !reverse );
}

/*! Return a filtered graph with only active edge (has at least one active contact)
 *
 * \tparam TGraph   Type of the graph (auto-deduced).
 * \param  graph    The graph to filter.
 * \return          A filtered graph (only the edges are filtered).
 */
template <
    typename TGraph
>
inline
filtered_graph <
    TGraph,
    std::function< bool( typename TGraph::edge_descriptor ) >
>
active_contact_graph( TGraph const& graph )
{
    return {
        graph,
        [&graph] ( typename TGraph::edge_descriptor e ) { 
            return std::any_of( graph[e].cbegin(), graph[e].cend(), is_active<decltype(graph[e][0])> );
        }
    };
}


/*! Return collision subgraphs with obstacle separation 
 *
 * It returns two separate component if they are only connected by an obstacle.
 *
 * \tparam TGraph   Type of the graph (auto-deduced).
 * \param  graph    The graph.
 * \return          A vector of subgraphs (each one represents a connected component).
 */
template < typename TGraph >
std::vector<TGraph>
collision_subgraphs( TGraph const& graph )
{
    // Obstacle ids
    std::vector<std::size_t> obstacle_ids;
    for ( auto const v : make_iterator_range( vertices(graph) ) )
        if ( graph[v]->is_obstacle() )
            obstacle_ids.push_back(v);

    // Filtered graph without obstacles
    auto const fgraph = graph_from_ids( graph, obstacle_ids, true );

    // Connected components of this filtered graph
    std::vector<std::size_t> components( num_vertices(graph) ); //! \warning clang 3.5.0 induced a bug while freeing this vector if constructed with {n} syntax. 
    const std::size_t num = connected_components( fgraph, &components[0] );

    // Eliminating obstacles from this list
    for ( std::size_t i : obstacle_ids )
        components[i] = num;

    // Initialization of the collision list
    std::vector<TGraph> subgraphs;
    subgraphs.reserve(num);

    // Foreach component (can be optimized ...)
    for ( std::size_t comp_id = 0; comp_id < num; ++comp_id )
    {
        // Ids of floes within this component
        std::vector<std::size_t> ids;
        for ( std::size_t i = 0; i < size(components); ++i )
            if ( components[i] == comp_id )
                ids.push_back(i);

        // Checking if some obstacles are connected with this component
        for ( std::size_t obs_id : obstacle_ids )
        {
            for ( auto const& v : make_iterator_range( adjacent_vertices(obs_id, graph) ) )
            {
                if ( components[v] == comp_id )
                {
                    ids.push_back(obs_id);
                    break;
                }
            }
        }

        // Creating sub graph if the component contains more than 1 floe|obstacle
        if ( ids.size() >= 2 )
        {
            subgraphs.push_back({});
            copy_graph( graph_from_ids( graph, ids ), subgraphs[subgraphs.size()-1]);
        }
    }
    
    return subgraphs;

}

/*! Return components with active contacts.
 *
 * \tparam TGraph   Type of the graph (auto-deduced).
 * \param graph     The contact graph.
 * \return          A vector of subgraphs extracted from the active contacts.
 */
template < typename TGraph >
std::vector<TGraph>
active_subgraphs( TGraph const& graph )
{
    // Filtered graph with only active edges
    const auto active_graph = active_contact_graph( graph );

    // Connected components of this filtered graph
    std::vector<std::size_t> components( num_vertices(graph) ); //! \warning clang 3.5.0 induced a bug while freeing this vector if constructed with {n} syntax. 
    const std::size_t num = connected_components( active_graph, &components[0] );

    // Initialization of the collision list
    std::vector<TGraph> subgraphs;
    subgraphs.reserve(num);

    // Foreach component (can be optimized ...)
    for ( std::size_t comp_id = 0; comp_id < num; ++comp_id )
    {
        // Ids of floes within this component
        std::vector<std::size_t> ids;
        for ( std::size_t i = 0; i < size(components); ++i )
            if ( components[i] == comp_id )
                ids.push_back(i);

        // Creating sub graph if the component contains more than 1 floe|obstacle
        if ( ids.size() >= 2 )
        {
            subgraphs.push_back({});

            // This version keeps all contacts between floes that have at least one active contact (even if it is with other floe)
            copy_graph( graph_from_ids( graph, ids ), subgraphs[subgraphs.size()-1]);

            // This version only keeps all contacts between two floes that have at least one active contact in common.
            //copy_graph( graph_from_ids( active_graph, ids ), subgraphs[subgraphs.size()-1]);
        }
    }
    
    return subgraphs;
}

/*! Return the number of contact in a graph
 *
 * \warning         The complexity is linear in the number of edges.
 * \tparam TGraph   Type of the graph (auto-deduced).
 * \param  graph    The graph.
 * \return          The number of contacts.
 */
template < typename TGraph >
std::size_t
num_contacts( TGraph const& graph )
{
    std::size_t cnt = 0;
    for ( auto const& e : make_iterator_range( edges(graph) ) )
        cnt += graph[e].size();
        
    return cnt;
}


/*! Return the number of contact in a graph
 *
 * \warning         The complexity is linear in the number of edges.
 * \tparam TGraph   Type of the graph (auto-deduced).
 * \param  graph    The graph.
 * \param  solved   Has the contact graph been solved ?
 * \return          The number of contacts.
 */
template < typename TGraph >
void
mark_solved( TGraph& graph, bool solved )
{
    for ( auto const& e : make_iterator_range( edges(graph) ) )
        graph[e].mark_solved(solved);
        
}

}} // namespace floe::collision

#endif // FLOE_COLLISION_CONTACT_GRAPH_HPP

