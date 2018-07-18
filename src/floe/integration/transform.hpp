/*!
 * \file floe/integration/transform.hpp
 * \brief Quadrature adapted to arbitrary shapes
 * \author Roland Denis
 * \date November 2014
 */

/*!
 * \namespace floe::integration::detail::transform
 * \brief Implementation details about transforming integration quadrature for a geometry.
 */

#ifndef FLOE_INTEGRATION_TRANSFORM_HPP_INCLUDED
#define FLOE_INTEGRATION_TRANSFORM_HPP_INCLUDED

#include <type_traits>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/int.hpp>

#include "floe/geometry/core/tag.hpp"
#include "floe/geometry/core/access.hpp"
#include "floe/geometry/core/coordinate_type.hpp"
#include "floe/geometry/core/static_num_points.hpp"
#include "floe/geometry/geometries/concepts/check.hpp"

#include "floe/integration/detail/multi_sum.hpp"


namespace floe { namespace integration {

namespace detail { namespace transform {
 
    /*! Type traits to get transform return type
     *
     * \tparam TFunction type of the function
     * \tparam TGeometry type ot the geometry
     */
    template <
        typename TFunction,
        typename TGeometry
    >
    struct return_type
    {
        //! The return type of a TFunction applied to a point of TGeometry
        typedef
            typename std::result_of< 
                TFunction&( 
                    typename floe::geometry::coordinate_type<TGeometry>::type,
                    typename floe::geometry::coordinate_type<TGeometry>::type
                )
            >::type
            type;
    };

    /*! Jacobian determinant for triangle transformation
     *
     * \tparam Triangle Triangle type
     * \param shape The triangle to which make the transformation
     */
    template < class Triangle >
    constexpr
    typename geometry::coordinate_type<Triangle>::type
    triangle_detJ ( const Triangle & shape )
    {
        using floe::geometry::get;
        return
              (get<0,0>(shape) - get<2,0>(shape)) * (get<1,1>(shape) - get<2,1>(shape))
            - (get<0,1>(shape) - get<2,1>(shape)) * (get<1,0>(shape) - get<2,0>(shape));
    }

    //! Quadrature transformation for a simple static polygon
    template <std::size_t N>
    struct simple_static_polygon_transform
    {
        BOOST_MPL_ASSERT_MSG
        (
            false, NOT_IMPLEMENTED_FOR_THIS_POLYGON, (types<boost::mpl::int_<N> >)
        );
    };
    
    //! Quadrature transformation for a triangle
    template <>
    struct simple_static_polygon_transform<3>
    {
        /*! Integration over an arbitrary triangle
         * \tparam TFunction Type of the function to integrate
         * \tparam TTriange  Type of the triangle on which to integrate
         * \tparam TStrategy Type of stragegy (quadrature) to be used
         * \param  f        The function to integrate
         * \param  shape    The triangle on which to integrate
         */
        template < 
            typename TFunction,
            typename TTriangle,
            typename TStrategy
        >
        static inline
        typename return_type<TFunction,TTriangle>::type
            apply( TFunction const& f, TTriangle const& triangle, TStrategy const& strategy ) 
        {
            using floe::geometry::get;
            //using T = typename return_type<TFunction,TTriangle>::type;
            using T = typename floe::geometry::coordinate_type<TTriangle>::type; 

            return triangle_detJ(triangle) * strategy.apply(
                [&f, &triangle] ( T x, T y ) {
                    return f(
                        x*get<0,0>(triangle) + y*get<1,0>(triangle) + (1-x-y)*get<2,0>(triangle),
                        x*get<0,1>(triangle) + y*get<1,1>(triangle) + (1-x-y)*get<2,1>(triangle)
                    );
                }
            );
        }

    };

}} // namespace detail::transform

namespace dispatch
{

//! Fallback transformation for unsupported geometries
template <
    typename TGeometry,
    typename Tag = typename floe::geometry::tag<TGeometry>::type
>
struct transform {
    BOOST_MPL_ASSERT_MSG
    (
        false, NOT_IMPLEMENTED_FOR_THIS_GEOMETRY_TYPE, (types<TGeometry>)
    );
};

/*! Quadrature transformation for a triangle
 *
 * Simply inherit of a SimpleStaticPolygon with 3 points
 */
template < typename TGeometry >
struct transform<TGeometry, floe::geometry::triangle_tag>
    : detail::transform::simple_static_polygon_transform<3>
{};

//! Quadrature transformation for simple static polygon
template < typename TGeometry >
struct transform<TGeometry, floe::geometry::simple_static_polygon_tag>
    : detail::transform::simple_static_polygon_transform< floe::geometry::static_num_points<TGeometry>::value >
{};

/*! Quadrature transformation for multiple simple static polygon
 *
 * It adds the result of the strategy applied on each SimpleStaticPolygon
 */
template < typename TMultiGeometry >
struct transform<TMultiGeometry, floe::geometry::multi_simple_static_polygon_tag>
    : detail::multi_sum
{
    template <
        typename TFunction,
        typename TStrategy
    >
    static inline 
    typename detail::transform::return_type<TFunction,TMultiGeometry>::type
        apply( TFunction const& function, TMultiGeometry const& multi, TStrategy const& strategy )
    {
        return multi_sum::apply
            <
                typename detail::transform::return_type<TFunction,TMultiGeometry>::type,
                transform<typename boost::range_value<TMultiGeometry>::type>
            >( function, multi, strategy );
    }
};

/*! Quadrature transformation for mesh
 *
 * It adds the result of the strategy applied on each cell of the mesh
 */
template < typename TMesh >
struct transform< TMesh, floe::geometry::mesh_tag >
{
    template <
        typename TFunction,
        typename TStrategy
    >
    static inline
    typename detail::transform::return_type<TFunction,TMesh>::type
        apply( TFunction const& function, TMesh const& mesh, TStrategy const& strategy )
    {
        return transform< typename floe::geometry::cells_type<TMesh>::type >::apply(
            function, floe::geometry::cells(mesh), strategy
        );
    }
};

} // namespace dispatch



/*! Transformation of a quadrature for an arbitrary shape.
 *
 * Transformation of a quadrature defined on a reference shape to be applied 
 * on an arbitrary shape.
 *
 * \tparam TFunction Type of the function to be integrate
 * \tparam TGeometry Type of the geometry
 * \tparam TStrategy Type of the quadrature method
 */
template <
    typename TFunction,
    typename TGeometry,
    typename TStrategy
>
inline
typename detail::transform::return_type<TFunction,TGeometry>::type
    transform( TFunction const& function, TGeometry const& geometry, TStrategy const& strategy )
{
    floe::geometry::concepts::check<TGeometry const>();
    return dispatch::transform<TGeometry>::apply( function, geometry, strategy );
}


}} // namespace floe::integration

#endif
