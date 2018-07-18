#ifndef FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_TRIANGLE_CONCEPT_HPP
#define FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_TRIANGLE_CONCEPT_HPP

#include <boost/concept_check.hpp>

#include "floe/geometry/geometries/concepts/simple_static_polygon_concept.hpp"
#include "floe/geometry/core/static_num_points.hpp"

namespace boost { namespace geometry { namespace concepts {

/*! Check triangle concept
 *
 * A triangle respect the SimpleStaticPolygon concept and has only 3 points.
 */
template < typename Geometry >
struct Triangle
{
    BOOST_CONCEPT_ASSERT( (concepts::SimpleStaticPolygon<Geometry>) );

public:

    BOOST_CONCEPT_USAGE(Triangle)
    {
        static_assert( geometry::static_num_points<Geometry>::value == 3, "A triangle must be made of 3 points." );
    }
};

//! Const version
template < typename Geometry >
struct ConstTriangle
{
    BOOST_CONCEPT_ASSERT( (concepts::ConstSimpleStaticPolygon<Geometry>) );

public:

    BOOST_CONCEPT_USAGE(ConstTriangle)
    {
        static_assert( geometry::static_num_points<Geometry>() == 3, "A triangle must be made of 3 points." );
    }
};


}}} // namespace boost::geometry::concepts

namespace floe { namespace geometry { namespace concepts
{
    using boost::geometry::concepts::Triangle;
    using boost::geometry::concepts::ConstTriangle;

}}} // namespace floe::geometry::concepts


#endif // FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_TRIANGLE_CONCEPT_HPP

