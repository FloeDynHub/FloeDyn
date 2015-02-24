/*!
 * \file floe/lcp/lcp.hpp
 * \brief Definition and manipulation of a LCP 
 * \author Roland Denis
 */

#ifndef FLOE_LCP_LCP_HPP
#define FLOE_LCP_LCP_HPP

#include <cstddef>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace floe { namespace lcp
{

/*! A LCP defined by A, q, w and z
 *
 * The problem is to find z that respect:
 * w = Az + q >= 0
 * z >= 0
 * w . z = 0
 *
 * \tparam T Fundamental type.
 */
template <
    typename T
>
struct LCP
{
    typedef boost::numeric::ublas::matrix<T> array_type;
    typedef boost::numeric::ublas::vector<T> vector_type;

    std::size_t dim;
    array_type A;
    vector_type q, w, z;

    LCP( std::size_t n )
        : dim(n), A(n, n), q(n), w(n), z(n)
    {}
};

/*! Return the solution error for a LCP
 *
 * \tparam T    Fundamental type (auto-deduced).
 * \param  lcp   The linear complementary problem.
 * \param  calc_w   True for recalculate w = Az + q
 */
template < typename T>
T LCP_error( LCP<T> const& lcp, bool calc_w = true )
{
    typedef typename LCP<T>::vector_type vector_type;
    using namespace boost::numeric::ublas;

    // Calculating w = Az + q
    vector_type w;
    if ( calc_w )
        w = prod(lcp.A, lcp.z) + lcp.q;
    else
        w = lcp.w;

    // Error on w
    T w_err = 0;
    for ( T value : w ) if (value < 0) w_err -= value;

    // Error on z
    T z_err = 0;
    for ( T value : lcp.z ) if (value < 0) z_err -= value;

    // Error on z.w
    T zw_err = 0;
    auto itz = lcp.z.begin();
    auto itw = w.begin();
    for ( ; itz != lcp.z.end(); ++itz, ++itw )
    {
        if (*itz != 0 && *itw != 0) 
            zw_err += std::abs( (*itz) * (*itw) );
    }

    return w_err + z_err + zw_err;
}

}} // namespace floe::lcp

#endif // FLOE_LCP_LCP_HPP

