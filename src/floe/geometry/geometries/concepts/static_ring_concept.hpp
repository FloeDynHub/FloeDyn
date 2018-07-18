#ifndef FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_STATIC_RING_CONCEPT_HPP
#define FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_STATIC_RING_CONCEPT_HPP

#include <type_traits>

#include <boost/concept_check.hpp>
#include <boost/range/concepts.hpp>
#include <boost/geometry/core/point_type.hpp>
#include <boost/geometry/core/coordinate_dimension.hpp>

#include "floe/geometry/core/access.hpp"

#include <boost/geometry/geometries/concepts/point_concept.hpp>

namespace boost { namespace geometry { namespace concept
{

/*! Static ring concept
 * 
 * It is like a ring but with compile-time fixed size.
 * The main advantage is that points can be accessed with get and set.
 * Therefore, compile-time expression can be used with such objects.
 */
template <typename Geometry>
class StaticRing
{
    typedef typename point_type<Geometry>::type point_type;

    BOOST_CONCEPT_ASSERT( (concept::Point<point_type>) );
    BOOST_CONCEPT_ASSERT( (boost::RandomAccessRangeConcept<Geometry>) );

    enum { pcount = dimension<Geometry>::value };

    template <
        typename R,
        std::size_t Dimension,
        std::size_t DimensionCount
    >
    struct dimension_checker
    {
        static void apply()
        {
            R* ring = 0;
            geometry::set<Dimension>(*ring, geometry::get<Dimension>(*ring));
            dimension_checker<R, Dimension+1, DimensionCount>::apply();
        }
    };

    template <
        typename R,
        std::size_t DimensionCount
    >
    struct dimension_checker<R, DimensionCount, DimensionCount>
    {
        static void apply() {}
    };

public:

    BOOST_CONCEPT_USAGE(StaticRing)
    {
        dimension_checker<Geometry, 0, pcount>::apply();
    }
};

template <typename Geometry>
class ConstStaticRing
{
    typedef typename point_type<Geometry>::type point_type;

    BOOST_CONCEPT_ASSERT( (concept::ConstPoint<point_type>) );
    BOOST_CONCEPT_ASSERT( (boost::RandomAccessRangeConcept<Geometry>) );

    enum { pcount = dimension<Geometry>::value };

    template <
        typename R,
        std::size_t Dimension,
        std::size_t DimensionCount
    >
    struct dimension_checker
    {
        static void apply()
        {
            const R* ring = 0;
            point_type point(geometry::get<Dimension>(*ring));
            boost::ignore_unused_variable_warning(point);
            dimension_checker<R, Dimension+1, DimensionCount>::apply();
        }
    };

    template <
        typename R,
        std::size_t DimensionCount
    >
    struct dimension_checker<R, DimensionCount, DimensionCount>
    {
        static void apply() {}
    };

public:

    BOOST_CONCEPT_USAGE(ConstStaticRing)
    {
        dimension_checker<Geometry, 0, pcount>::apply();
    }
};

}}} // namespace boost::geometry::concepts


#endif // FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_STATIC_RING_CONCEPT_HPP

