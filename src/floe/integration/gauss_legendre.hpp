/*!
 * \file gauss_legendre.hpp 
 * \brief Gauss-Legendre integration methods
 * \author Roland Denis
 */

#ifndef INTEGRATION_GAUSS_LEGENDRE_HPP_INCLUDED
#define INTEGRATION_GAUSS_LEGENDRE_HPP_INCLUDED

#include <type_traits>

#include "floe/integration/quadrature.hpp"

namespace floe { namespace integration {

/*! Gauss-Legendre integration method
 *
 * \tparam T Space fondamental type
 * \tparam D Space dimension
 * \tparam N Method order
 *
 * Gauss-Legendre integration method archetype 
 * useable on reference shape (eg unit interval, triangle, ...)
 */
template <
    typename T,     // Space fondamental type
    unsigned int D, // Space dimension
    unsigned int N  // Method order
>
class RefGaussLegendre;

/*! Gauss-Legendre integration method, order 2, 2D space
 *
 * Implementation of order 2 on 2D space with the reference triangle defined by (0,0), (1,0) and (1,1).
 *
 * \tparam T Space fondamental type
 *
 * \todo   Verify order
 *
 * \remark Why only triangle ? Why not other shapes ?
 */
template <
    typename T  // Space fondamental type
>
class RefGaussLegendre<T, 2, 2> 
    : public Quadrature< RefGaussLegendre<T,2,2>, T, 2 >
{
    //! Integration points in clockwise order
    static const T p1[2]; 
    static const T p2[2];
    static const T p3[2];

    //! Integration weights
    static const T weight[3];

    
public :
    
    /*! Integration of a function on the reference triangle
     *
     * \tparam Function Type of the function to integrate
     * \param f Function to integrate
     * \result The integration of f on the reference triangle
     *
     * \remark Method that can be mutualized between implementations.
     * It is why array are used instead of inlined values
     */
    template < typename Function >
    constexpr
    auto apply_impl( const Function & f ) const
        -> typename std::result_of<Function&(T,T)>::type
    {
        return 
              weight[0] * f( p1[0], p1[1] )
            + weight[1] * f( p2[0], p2[1] )
            + weight[2] * f( p3[0], p3[1] )
        ;
    }

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> pointsAndWeights()  
    {
        Eigen::Matrix<T, 3, 3> m;
        m << p1[0], p1[1], weight[0],
            p2[0], p2[1], weight[1],
            p3[0], p3[1], weight[2];
        return m;
    }
};

template <typename T>
const T RefGaussLegendre<T,2,2>::p1[2] = { T(1)/T(6), T(1)/T(6) };

template <typename T>
const T RefGaussLegendre<T,2,2>::p2[2] = { T(2)/T(3), T(1)/T(6) };

template <typename T>
const T RefGaussLegendre<T,2,2>::p3[2] = { T(1)/T(6), T(2)/T(3) };

template <typename T>
const T RefGaussLegendre<T,2,2>::weight[3] = { T(1)/T(6), T(1)/T(6), T(1)/T(6) };


/*! Gauss-Legendre integration method, order 1, 2D space
 *
 * Implementation of order 1 on 2D space with the reference triangle defined by (0,0), (1,0) and (1,1).
 *
 * \tparam T Space fondamental type
 *
 * \todo   Verify order
 *
 * \remark Why only triangle ? Why not other shapes ?
 */
template <
    typename T  // Space fondamental type
>
class RefGaussLegendre<T, 2, 1> 
    : public Quadrature< RefGaussLegendre<T,2,1>, T, 2 >
{
    //! Integration points in clockwise order
    static const T p1[2]; 

    //! Integration weights
    static const T weight;

    
public :
    
    /*! Integration of a function on the reference triangle
     *
     * \tparam Function Type of the function to integrate
     * \param f Function to integrate
     * \result The integration of f on the reference triangle
     *
     * \remark Method that can be mutualized between implementations.
     * It is why array are used instead of inlined values
     */
    template < typename Function >
    constexpr
    auto apply_impl( const Function & f ) const
        -> typename std::result_of<Function&(T,T)>::type
    {
        return weight * f( p1[0], p1[1] );
    }
};

template <typename T>
const T RefGaussLegendre<T,2,1>::p1[2] = { T(1)/T(3), T(1)/T(3) };

template <typename T>
const T RefGaussLegendre<T,2,1>::weight = T(1)/T(2);

}} // namespace floe::integration

#endif
