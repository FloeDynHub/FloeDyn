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
    typedef boost::numeric::ublas::matrix<T> array_type;    //!< Type of array.
    typedef boost::numeric::ublas::vector<T> vector_type;   //!< Type of vector.

    std::size_t dim;    //!< Dimension of the problem. (4 times the number of contacts)
    array_type A;       //!< The matrix of the problem.
    vector_type q;      //!< The vector of the problem.
    vector_type w;      //!< The vector complementary to the LCP solution.
    vector_type z;      //!< The solution of the LCP.

    //! Constructor given the dimension of the problem.
    LCP( std::size_t n ) : dim(n), A(n, n), q(n), w(n), z(n,0.0) {}
};

/*! Return the solution error for a LCP
 *
 * \tparam T    Fundamental type (auto-deduced).
 * \param  lcp   The linear complementary problem.
 * \param  calc_w   True for recalculate w = Az + q
 * \return |w^-|_0 + |z^-|_0 + |w.z|
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

template < typename T>
boost::numeric::ublas::vector<T> LCP_error_detailed( LCP<T> const& lcp)
{
    typedef typename LCP<T>::vector_type vector_type;
    using namespace boost::numeric::ublas;

    std::size_t dim = lcp.dim;
    vector_type Vec_Err(3*dim,0);

    // Calculation: w = Az + q;
    vector_type w;
    w = prod(lcp.A, lcp.z) + lcp.q;

    // Calculation: zw = z^T w;
    vector_type zw(dim);
    for (std::size_t i=0; i<dim; ++i) {
        zw(i) = lcp.z(i) * w(i);
    }

    for (std::size_t i=0; i<dim; ++i) {
        Vec_Err(i) = zw(i); // energy part 
        if (lcp.z(i)<0) {Vec_Err(i+dim) = -lcp.z(i);} // impulse part
        if (w(i)<0) {Vec_Err(i+2*dim) = -w(i);} // relative velocities after contact
    }

    return Vec_Err;
}

template < typename T>
T LCP_error_global( boost::numeric::ublas::vector<T> const& Err )
{
    T LCP_err(0);
    // global LCP error:
    for (std::size_t i=0; i<Err.size(); ++i) {
        LCP_err += std::abs(Err(i));
    }

    return LCP_err;
}

// Calc error new to highlight the coulomb failure:
// ajout Matthias
template < typename T>
boost::numeric::ublas::vector<T> LCP_err_Coul( LCP<T> const& lcp)
{
    typedef typename LCP<T>::vector_type vector_type;

    int m;
    m = lcp.dim/4;
    // std::cout << "m: " << m << ",\n";

    using namespace boost::numeric::ublas;
    compressed_matrix<T, column_major> E;
    // Friction coupling matrix E
    E.resize(2*m, m, false);
    for (int j = 0; j < m; ++j ) {
        E(2*j+1, j) = E(2*j, j) = 1;
    }
    E = trans(E);

    vector_type cond_Coulomb(m);
    for (int i=0; i<m; ++i) {
        int j=2*i; int k=2*i+1;
        // std::cout   << "i= " << i << ",\n z(i): "
        //             << lcp.z(i) << ",\n E(i,j): "
        //             << E(i,j) << ", E(i,k): " << E(i,k) << ",\n z(j+m): "
        //             << lcp.z(j+m) << ", z(k+m): " << lcp.z(k+m) << ",\n" << std::endl;
        cond_Coulomb(i) = 0.7*lcp.z(i)-E(i,j)*lcp.z(j+m)-E(i,k)*lcp.z(k+m);
    }

    return cond_Coulomb;
}
//endmatt

}} // namespace floe::lcp

#endif // FLOE_LCP_LCP_HPP

