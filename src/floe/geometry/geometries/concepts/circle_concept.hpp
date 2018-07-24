#ifndef FLOE_GEOMETRY_GEOMETRIES_CONCEPT_CIRCLE_CONCEPT_HPP
#define FLOE_GEOMETRY_GEOMETRIES_CONCEPT_CIRCLE_CONCEPT_HPP

#include <cstddef>

#include <boost/concept_check.hpp>

#include <boost/geometry/core/point_type.hpp>

#include "floe/geometry/core/access.hpp"
#include "floe/geometry/core/radius.hpp"
#include "floe/geometry/core/coordinate_dimension.hpp"

namespace boost { namespace geometry { namespace concepts {

/*! Checks circle concept (const version)
 */
template < typename Geometry >
class ConstCircle
{
    typedef typename geometry::point_type<Geometry>::type  point_type;
    typedef typename geometry::radius_type<Geometry>::type radius_type;

    template < std::size_t Dimension, std::size_t DimensionCount >
    struct dimension_checker
    {
        static void apply()
        {
            typedef typename coordinate_type<Geometry>::type coordinate_type;
            const Geometry* s = 0;
            coordinate_type coord(geometry::get<Dimension>(*s));
            boost::ignore_unused_variable_warning(coord);
            dimension_checker<Dimension + 1, DimensionCount>::apply();
        }
    };

    template <std::size_t DimensionCount>
    struct dimension_checker<DimensionCount, DimensionCount>
    {
        static void apply() {}
    };

public:

    BOOST_CONCEPT_USAGE(ConstCircle)
    {
        static const std::size_t n = dimension<Geometry>::value;
        static_assert( n == 2, "circles can only be of dimension 2" );
        dimension_checker<0, 2>::apply();
        dimension_checker<0, 2>::apply(); // Why 2 times ?

        // Check radius access
        Geometry const* s = 0;
        radius_type coord(floe::geometry::get_radius(*s));
        boost::ignore_unused_variable_warning(coord);
    }
};


/*! Checks circle concept (non const version)
 */
template < typename Geometry >
class Circle
{
    BOOST_CONCEPT_ASSERT( (concepts::ConstCircle<Geometry>) );

    typedef typename geometry::point_type<Geometry>::type point_type;
    typedef typename geometry::radius_type<Geometry>::type radius_type;

    template < std::size_t Dimension, std::size_t DimensionCount >
    struct dimension_checker
    {
        static void apply()
        {
            Geometry* s = 0;
            geometry::set<Dimension>(*s, geometry::get<Dimension>(*s) );
            dimension_checker<Dimension + 1, DimensionCount>::apply();
        }
    };

    template <std::size_t DimensionCount>
    struct dimension_checker<DimensionCount, DimensionCount>
    {
        static void apply() {}
    };

public:

    BOOST_CONCEPT_USAGE(Circle)
    {
        static const std::size_t n = dimension<Geometry>::value;
        static_assert( n == 2, "circles can only be of dimension 2");
        dimension_checker<0, 2>::apply();
        dimension_checker<0, 2>::apply(); // Why 2 times ?

        // Check radius access
        Geometry* s = 0;
        set_radius(*s, get_radius(*s));
    }
};

}}} // namespace boost::geometry::concepts

namespace floe { namespace geometry { namespace concepts
{
    using boost::geometry::concepts::Circle;
    using boost::geometry::concepts::ConstCircle;

}}} // namespace floe::geometry::concepts

#endif // FLOE_GEOMETRY_GEOMETRIES_CONCEPT_CIRCLE_CONCEPT_HPP
