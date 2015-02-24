/*!
 * \file floe/geometry/algorithms/transform.hpp
 * \brief Transform (translate, rotate, scale, ...) for circles and triangle mesh. Also imports boost versions.
 * \author Roland Denis
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/geometry/reference/algorithms/transform.html">boost/geometry/algorithms/transform.hpp</a>
 *
 * \namespace boost::geometry::detail::transform
 * \brief Implementation details for transform calculation.
 */

#ifndef FLOE_GEOMETRY_ALGORITHMS_TRANSFORM
#define FLOE_GEOMETRY_ALGORITHMS_TRANSFORM

#include <boost/geometry/algorithms/transform.hpp>
#include <boost/geometry/algorithms/assign.hpp>

#include "floe/geometry/geometries/triangle_mesh.hpp"
#include "floe/geometry/core/radius.hpp"
#include "floe/geometry/views/center_view.hpp"

namespace boost { namespace geometry
{

namespace detail { namespace transform
{

namespace fg = floe::geometry;

/*! Transformation with strategy for TriangleMesh
 *
 * There is a specialization only for floe::geometry::TriangleMesh
 * \todo a better solution ...
 */
struct transform_TriangleMesh
{

    template <
        typename TPoint1,
        typename TPoint2,
        typename Strategy
    >
    static inline 
    bool apply(
        TriangleMesh<TPoint1> const& triangle_mesh1,
        TriangleMesh<TPoint2> & triangle_mesh2,
        Strategy const& strategy
    )
    {
        // Copy of the connectivity
        triangle_mesh2.connectivity() = triangle_mesh1.connectivity();
        
        // Transformation of the points
        return boost::geometry::transform( triangle_mesh1.points(), triangle_mesh2.points(), strategy );
    }
};

/*! Transformation with strategy for circle_tag geometry
 *
 * \warning The radius is not scaled !!!
 */
struct transform_circle
{
    template <
        typename TCircle1,
        typename TCircle2,
        typename TStrategy
    >
    static inline
    bool apply( TCircle1 const& circle1, TCircle2& circle2, TStrategy const& strategy )
    {
        set_radius(circle2, get_radius(circle1));
        
        center_view<TCircle2> center2(circle2);
        return boost::geometry::transform(
            center_view<const TCircle1>(circle1),
            center2,
            strategy
        );
    }

};

}} // namespace detail::transform

namespace dispatch
{
    template <
        typename TPoint1,
        typename TPoint2
    >
    struct transform< TriangleMesh<TPoint1>, TriangleMesh<TPoint2>, mesh_tag, mesh_tag >
        : detail::transform::transform_TriangleMesh
    {};

    template <
        typename TCircle1,
        typename TCircle2
    >
    struct transform< TCircle1, TCircle2, circle_tag, circle_tag >
        : detail::transform::transform_circle
    {};

} // namespace dispatch

}} // namespace boost::geometry

namespace floe { namespace geometry
{
    using boost::geometry::transform;

}} // namespace floe::geometry


#endif // FLOE_GEOMETRY_ALGORITHMS_TRANSFORM
