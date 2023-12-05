/*!
 * \file floe/collision/contact_graph.hpp
 * \brief Contact graph definition and manipulation.
 * \author Roland Denis, Quentin Jouet
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
#include "floe/collision/floe_vertex.hpp"

#include <iostream> // DEBUG
// #include <chrono> // test

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
    // typename TContactPoint::floe_type const*, // Vertex Property
    FloeVertex<typename TContactPoint::floe_type>, // Vertex Property
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
        [&graph] ( typename TGraph::vertex_descriptor v ) { return ! graph[v].floe->is_obstacle(); } // fv_test
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


template <
    typename TGraph
>
inline
filtered_graph <
    TGraph,
    boost::keep_all,
    std::function< bool( typename TGraph::vertex_descriptor ) >
>
no_lone_node_subgraph( TGraph const& graph )
{
    return {
        graph,
        boost::keep_all(),
        [&graph] ( typename TGraph::vertex_descriptor v ) {
            return (in_degree(v, graph) != 0 || out_degree(v, graph) != 0);
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
        if ( graph[v].floe->is_obstacle() ) // fv_test
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
        for ( std::size_t i = 0; i < components.size(); ++i )
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
    const auto active_graph1 = active_contact_graph( graph );
    std::vector<std::size_t> active_vertex_ids;
    for ( auto const& edge : make_iterator_range( edges( active_graph1 ) ) )
    {
        active_vertex_ids.push_back(source(edge, active_graph1));
        active_vertex_ids.push_back(target(edge, active_graph1));

    }
    if (active_vertex_ids.size() == 0) return {};
    auto const active_graph = graph_from_ids( active_graph1, active_vertex_ids );

    // Connected components of this filtered graph
    std::vector<int> components( num_vertices(graph), -1 ); //! \warning clang 3.5.0 induced a bug while freeing this vector if constructed with {n} syntax. 
    const std::size_t num = connected_components( active_graph, &components[0] );
    std::vector<std::array<int, 2>> active_vertex_map;
    for (int i = 0; std::size_t(i) < components.size(); ++i)
        if (components[i] != -1) active_vertex_map.push_back({{i, components[i]}});
    
    // Initialization of the collision list
    std::vector<TGraph> subgraphs;
    subgraphs.reserve(num);

    // Foreach component
    for ( std::size_t comp_id = 0; comp_id < num; ++comp_id )
    {
        // Ids of floes within this component
        std::vector<std::size_t> ids;
        for ( auto e : active_vertex_map)
            if ( std::size_t(e[1]) == comp_id )
                ids.push_back(e[0]);

            subgraphs.push_back({});
            auto& subgraph = subgraphs.back();

            // This version keeps all contacts between floes that have at least one active contact (even if it is with other floe)
            copy_graph( graph_from_ids( graph, ids ), subgraph);
            // auto nb_contact = num_contacts(subgraph);
            // if ( nb_contact > 100 ){
            //     subgraphs.pop_back();
            //     subgraphs.push_back({});
            //     auto& subgraph = subgraphs.back();
            // This version only keeps all contacts between two floes that have at least one active contact in common.
            //     copy_graph( graph_from_ids( active_graph, ids ), subgraph);
            //     std::cout << "Reduce graph " <<  nb_contact << " -> " << num_contacts(subgraph) << std::endl;
            // }

            // Keep track of ids in parent graph (for optimisations)
            std::size_t i = 0;
            for (auto const& v : make_iterator_range( vertices( subgraph ) )){
                subgraph[v].parent_descriptor = ids[i]; ++i;
            }

            // This version only keeps all contacts between two floes that have at least one active contact in common.
            //copy_graph( graph_from_ids( active_graph, ids ), subgraphs[subgraphs.size()-1]);
    }
    
    return subgraphs;
}


template < typename TGraph >
std::vector<TGraph>
quad_cut( TGraph const& graph )
{
    using real_type = double;
    // Calc window
    real_type mg = 1e-8; // margin
    real_type min_x, min_y, max_x, max_y;
    min_x = min_y = std::numeric_limits<real_type>::max();
    max_x = max_y = - std::numeric_limits<real_type>::max();
    for ( auto const v : make_iterator_range( vertices(graph) ) ){
        max_x = std::max(graph[v].floe->state().pos.x, max_x);
        min_x = std::min(graph[v].floe->state().pos.x, min_x);
        max_y = std::max(graph[v].floe->state().pos.y, max_y);
        min_y = std::min(graph[v].floe->state().pos.y, min_y);
    }
    std::array<real_type, 4> win = {{ min_x - mg, max_x + mg, min_y - mg, max_y + mg }};
    // make ids vectors
    int side_len = 2;
    std::vector<std::vector<std::size_t>> idss(side_len * side_len);
    real_type parcel_width = (win[1] - win[0]) / side_len;
    real_type parcel_height = (win[3] - win[2]) / side_len;
    for ( auto const v : make_iterator_range( vertices(graph) ) ){
        auto XX = (graph[v].floe->state().pos.x - win[0]) / parcel_width;
        auto YY = (graph[v].floe->state().pos.y - win[2]) / parcel_height;
        int XID{(int)floor(XX)};
        int YID{(int)floor(YY)};
        int sub_graph_id = YID * side_len + XID;
        idss[sub_graph_id].push_back(v);
    }

    std::vector<TGraph> subgraphs;
    // Foreach ids vector
    for ( auto& ids : idss )
    {
            subgraphs.push_back({});
            auto& subgraph = subgraphs.back();

            // This version keeps all contacts between floes that have at least one active contact (even if it is with other floe)
            copy_graph( graph_from_ids( graph, ids ), subgraph);

            // Keep track of ids in parent graph (for optimisations)
            std::size_t i = 0;
            for (auto const& v : make_iterator_range( vertices( subgraph ) )){
                subgraph[v].parent_descriptor = ids[i]; ++i;
            }
    }
    
    return subgraphs;
}


/*! Return the number of contact in a graph
 *
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


/*! Mark contacts as solved
 *
 * \tparam TGraph   Type of the graph (auto-deduced).
 * \param  graph    The graph.
 * \param  solved   Has the contact graph been solved ?
 * \return          void
 */
template < typename TGraph >
void
mark_solved( TGraph& graph, bool solved )
{
    for ( auto const& e : make_iterator_range( edges(graph) ) )
        graph[e].mark_solved(solved);
}


/*! Mark contacts as changed, i.e need relative_speed refresh (cache optimization)
 *
 * \tparam TGraph   Type of the graph (auto-deduced).
 * \param  graph    The graph.
 * \return          void
 */
template < typename TGraph >
void
mark_changed( TGraph& graph )
{
    for ( auto const& e : make_iterator_range( edges(graph) ) )
        graph[e].mark_changed();
}


/*! Mark contacts as changed in parent graph, i.e need relative_speed refresh (cache optimization)
 * for all contacts with floes (vertices out-edges) of the child graph
 *
 * \tparam TGraph   Type of the graph (auto-deduced).
 * \param  graph    The graph.
 * \return          void
 */
template < typename TGraph >
void
mark_changed_parent( TGraph& graph, TGraph& parent_graph )
{
    for ( auto const v : boost::make_iterator_range( vertices(graph) ) ){
        for (auto const e : make_iterator_range( out_edges(graph[v].parent_descriptor, parent_graph) ) )
            parent_graph[e].mark_changed();
    }
    
}


/*! Return the number of contact in a graph or a filtered graph
 *
 * \tparam TGraph   Type of the graph (auto-deduced).
 * \param  graph    The graph.
 * \return          The number of contacts.
 */
template < typename TGraph >
std::size_t 
filtered_num_vertices( TGraph& graph )
{
    std::size_t i = 0;
    for ( auto const& e : make_iterator_range( edges(graph) ) ) ++i;
    return i;
        
}

}} // namespace floe::collision

#endif // FLOE_COLLISION_CONTACT_GRAPH_HPP

