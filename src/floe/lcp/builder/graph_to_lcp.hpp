/*!
 * \file floe/lcp/builder/graph_to_lcp.hpp
 * \brief Contact graph to LCP converter.
 * \author Roland Denis
 */

#ifndef FLOE_LCP_BUILDER_GRAPH_TO_LCP_HPP
#define FLOE_LCP_BUILDER_GRAPH_TO_LCP_HPP

#include <cstddef>
#include <algorithm>
#include <utility>
#include <type_traits>

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation_blocked.hpp>
#include <boost/numeric/ublas/vector.hpp>
// #include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/graph/graph_utility.hpp>

#include "floe/geometry/geometry.hpp" // TODO lighter include
#include "floe/geometry/core/access.hpp"
#include "floe/geometry/arithmetic/arithmetic.hpp"
#include "floe/geometry/arithmetic/determinant.hpp"
#include "floe/lcp/lcp.hpp"

namespace floe { namespace lcp { namespace builder
{

namespace ublas = boost::numeric::ublas;


/*! LCP build from a contact graph.
 *
 * \tparam T        Fundamental type.
 * \tparam TGraph   Type of the graph.
 *
 * A reference to the graph is keeped to be able, later, to return a restriction of the LCP
 * when a subgraph is given.
 */
template <
    typename T,
    typename TGraph
>
class GraphLCP
{

public:

    //! Delete default constructor
    GraphLCP() = delete;

    //! Constructor from a graph
    GraphLCP( TGraph const& graph ) : m_graph(graph)
    {
        init();
    }

    //! Initialize internal matrices need to build a LCP.
    void init();

    //! Return the main matrix of the LCP
    //! \todo parameter for decompression
    lcp::LCP<T> getLCP() const;
    lcp::LCP<T> getLCP_d(lcp::LCP<T> const& lcp_c, ublas::vector<T> const& Solc, T epsilon = 0) const;
    T born_sup_d(lcp::LCP<T> const& lcp_c, T epsilon = 0) const;

    ublas::diagonal_matrix<T>   M;   //!< Mass and momentum matrix.
    ublas::diagonal_matrix<T>   invM;  //!< Inverse of the mass/momentum matrix.
    ublas::compressed_matrix<T, ublas::column_major> J;      //!< Normal matrix.
    ublas::compressed_matrix<T, ublas::column_major> D;      //!< Tangantial matrix.
    ublas::compressed_matrix<T, ublas::column_major> E;      //!< Friction coupling matrix.
    ublas::diagonal_matrix<T>   mu;     //<! Static friction coefficients.
    mutable ublas::vector<T> W; //<! Floes' speed vector.

private:
    TGraph const& m_graph; //<! The contact graph.
};

template < typename T, typename TGraph >
void
GraphLCP<T, TGraph>::
init()
{
    using namespace std;
    namespace fg = floe::geometry;
    
    // Number of floes 
    const size_t n = num_vertices( m_graph );

    // Number of contact
    const size_t m = num_contacts( m_graph );
    
    // Mass & momentum matrix (and inverse) initialization
    M.resize(3*n, 3*n, false);
    invM.resize(3*n, 3*n, false);

    for ( size_t i = 0; i < 3*n; i += 3 )
    {
        if ( m_graph[i/3]->is_obstacle() )
        {
            M(i+2,i+2)  = M(i+1, i+1)  = M(i, i) = std::numeric_limits<T>::max();
            invM(i+2,i+2) = invM(i+1, i+1) = invM(i, i) = 0;
        } 
        else 
        {
            M(i+1, i+1) = M(i, i) = m_graph[i/3]->mass();
            M(i+2, i+2) = m_graph[i/3]->moment_cst();
            invM(i+1, i+1) = invM(i, i) = T(1) / m_graph[i/3]->mass();
            invM(i+2, i+2) = T(1) / m_graph[i/3]->moment_cst();
        }
    }

            
    // Normal matrix J, tangential matrix D and static friction matrix
    J.resize(3*n, m, false);
    D.resize(3*n, 2*m, false);
    mu.resize(m, m, false);

    size_t j = 0;

    // Foreach edges of the graph ...
    for ( auto const& edge : make_iterator_range( edges( m_graph ) ) )
    {

        const pair<size_t, size_t> id = minmax( source(edge, m_graph), target(edge, m_graph) );
        auto const* floe1 = m_graph[id.first];
        auto const* floe2 = m_graph[id.second];

        const size_t i1 = 3*id.first;
        const size_t i2 = 3*id.second;

        // Foreach contact between this 2 floes ...
        for ( auto const& contact : m_graph[edge] )
        {
            const T c = ( floe1 == contact.floe1 ) ? -1 : 1;

            typedef typename std::decay<decltype(contact)>::type::point_type point_type;

            // Normal and tangent to the contact
            const point_type normal  = contact.frame.v();
            const point_type tangent = contact.frame.u();

            // Vector from the mass center to the contact point
            point_type r1 = contact.frame.center(); fg::subtract_point( r1, floe1->frame().center() );
            point_type r2 = contact.frame.center(); fg::subtract_point( r2, floe2->frame().center() );

            // Filling normal matrix
            J(i1, j)   = c * fg::get<0>(normal);
            J(i1+1, j) = c * fg::get<1>(normal);
            J(i1+2, j) = c * fg::determinant<T>( r1, normal );

            J(i2, j)   = -c * fg::get<0>(normal);
            J(i2+1, j) = -c * fg::get<1>(normal);
            J(i2+2, j) = -c * fg::determinant<T>( r2, normal );

            // Filling tangential matrix
            //      first part:
            D(i1, 2*j)   = c * fg::get<0>(tangent);
            D(i1+1, 2*j) = c * fg::get<1>(tangent);
            D(i1+2, 2*j) = c * fg::determinant<T>( r1, tangent );

            D(i2, 2*j)   = -c * fg::get<0>(tangent);
            D(i2+1, 2*j) = -c * fg::get<1>(tangent);
            D(i2+2, 2*j) = -c * fg::determinant<T>( r2, tangent );

            //      second part:
            D(i1, 2*j+1)   = -c * fg::get<0>(tangent);
            D(i1+1, 2*j+1) = -c * fg::get<1>(tangent);
            D(i1+2, 2*j+1) = -c * fg::determinant<T>( r1, tangent );

            D(i2, 2*j+1)   = c * fg::get<0>(tangent);
            D(i2+1, 2*j+1) = c * fg::get<1>(tangent);
            D(i2+2, 2*j+1) = c * fg::determinant<T>( r2, tangent );

            // Filling static friction coefficient matrix
            mu(j, j) = floe1->mu_static(); // An expression depending of floe1 & floe2 mu ?
            
            ++j;
        } // foreach contact.
    } // foreach edges.

    // Friction coupling matrix E
    E.resize(2*m, m, false);
    for ( j = 0; j < m; ++j )
        E(2*j+1, j) = E(2*j, j) = 1;

}

template <typename T, typename TGraph>
floe::lcp::LCP<T>
GraphLCP<T, TGraph>::
getLCP() const
{
    using namespace boost::numeric::ublas;
    
    const std::size_t n = J.size1()/3; // Hum ...
    const std::size_t m = J.size2();
    
    // LCP main matrix
    lcp::LCP<T> lcp(4*m);
    auto & A = lcp.A;

    // Temporaries
    decltype(J) iMJ(3*n, m);        axpy_prod(invM, J, iMJ, true);
    decltype(D) iMD(3*n, 2*m);    axpy_prod(invM, D, iMD, true);

    // Filling by blocks
    //      1st row
    
    //axpy_prod( trans(J), iMJ, project( A, range(0, m), range(0, m) ),   true );
    project(A, range(0, m), range(0, m)) = prod(trans(J), iMJ);
    //matrix<T> block11(m,m); axpy_prod( trans(J), iMJ, block11, true ); project(A, range(0,m), range(0, m)) = block11;
    //matrix<T> block11 = block_prod<matrix<T>,16>( trans(J), iMJ); project(A, range(0,m), range(0, m)) = block11;
    
    //axpy_prod( trans(J), iMD, project( A, range(0, m), range(m, 3*m) ), true );
    project(A, range(0, m), range(m, 3*m)) = prod(trans(J), iMD);
    //matrix<T> block12(m, 2*m); axpy_prod( trans(J), iMD, block12, true); project(A, range(0,m), range(m,3*m)) = block12;
    //matrix<T> block12 = block_prod<matrix<T>,16>( trans(J), iMD); project(A, range(0,m), range(m,3*m)) = block12;

    project( A, range(0, m), range(3*m, 4*m) ) = zero_matrix<T>(m, m);

    //      2nd row
    
    //axpy_prod( trans(D), iMJ, project( A, range(m, 3*m), range(0, m) ),   true );
    //project(A, range(m, 3*m), range(0, m)) = prod(trans(D), iMJ);
    project(A, range(m, 3*m), range(0, m)) = trans(project(A, range(0, m), range(m, 3*m))); // It is faster to transpose already calculated block
    //project(A, range(m, 3*m), range(0, m)) = trans(block12);
    
    //axpy_prod( trand(D), iMD, project( A, range(m, 3*m), range(m, 3*m) ), true );
    project(A, range(m, 3*m), range(m, 3*m)) = prod(trans(D), iMD);
    //matrix<T> block22(2*m, 2*m); axpy_prod( trans(D), iMD, block22, true); project(A, range(m, 3*m), range(m, 3*m)) = block22;
    //matrix<T> block22 = block_prod<matrix<T>,16>( trans(D), iMD); project(A, range(m, 3*m), range(m, 3*m)) = block22;

    project( A, range(m, 3*m), range(3*m, 4*m) ) = E;

    //      3th row
    project( A, range(3*m, 4*m), range(0, m) ) = mu;
    project( A, range(3*m, 4*m), range(m, 3*m) ) = -trans(E);
    project( A, range(3*m, 4*m), range(3*m, 4*m) ) = zero_matrix<T>(m,m);


    // And now, the q vector !!
    // vector<T> W(3*n); // changed to access W from outside
    W.resize(3*n, false);

    // The speed vector is not prepared before to be sync with the actual floes states
    namespace fg = floe::geometry;
    for (std::size_t i = 0; i < 3*n; i+=3)
    {
        auto const state = m_graph[i/3]->state();
        W(i)   = fg::get<0>(state.speed);
        W(i+1) = fg::get<1>(state.speed);
        W(i+2) = state.rot;
    }

    // Filling q
    auto& q = lcp.q;
    project(q, range(0, m))   = prod( trans(J), W );
    project(q, range(m, 3*m)) = prod( trans(D), W );
    project(q, range(3*m, 4*m)) = zero_vector<T>(m);

    // Job done !
    return lcp;
}


template <typename T, typename TGraph>
floe::lcp::LCP<T>
GraphLCP<T, TGraph>::
getLCP_d(lcp::LCP<T> const& lcp_c, ublas::vector<T> const& Solc, T epsilon) const
{
    using namespace boost::numeric::ublas;
    const std::size_t m = J.size2();

    auto lcp_d = lcp_c;
    ublas::vector<T> ezc = epsilon * subrange(lcp_c.z, 0, m);

    // Filling q_d
    auto& q = lcp_d.q;
    
    ublas::vector<T> temp2 = prod( invM, ublas::vector<T>(prod(J, ezc)) );
    project(q, range(0, m)) = prod( trans(J), temp2 ) + prod( trans(J), Solc );
    project(q, range(m, 3*m)) = prod( trans(D), temp2 ) + prod( trans(D), Solc);
    project(q, range(3*m, 4*m)) = zero_vector<T>(m);

    return lcp_d;
}

template <typename T, typename TGraph>
T
GraphLCP<T, TGraph>::
born_sup_d(lcp::LCP<T> const& lcp_c, T epsilon) const
{
    using namespace boost::numeric::ublas;
    const std::size_t m = J.size2();

    ublas::vector<T> ezc = epsilon * subrange(lcp_c.z, 0, m);
    return inner_prod(
            ezc,
            prod(
                trans(J),
                ublas::vector<T>(prod(invM,
                    ublas::vector<T>(prod(J, ezc))
                ))
            ) / inner_prod( prod(trans(W), M), W) );
}


}}} // namespace floe::lcp::builder

#endif // FLOE_LCP_BUILDER_GRAPH_TO_LCP_HPP

