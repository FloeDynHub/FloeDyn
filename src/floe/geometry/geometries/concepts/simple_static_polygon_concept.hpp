#ifndef FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_SIMPLE_STATIC_POLYGON_CONCEPT_HPP
#define FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_SIMPLE_STATIC_POLYGON_CONCEPT_HPP

#include <cstddef>      // For std::size_t
#include <type_traits>  // std::decay

#include <boost/concept_check.hpp>
#include <boost/range/concepts.hpp>

#include <boost/geometry/core/point_type.hpp>

#include "floe/geometry/core/access.hpp"
#include "floe/geometry/core/static_num_points.hpp"

#include <boost/geometry/geometries/concepts/point_concept.hpp>

namespace boost { namespace geometry { namespace concept {

/*! Checks simple static polygon concept
 *
 * A simple static polygon is a polygon made of a fixed
 * number of points and composed only by an non-closed exterior ring.
 * It can be, for example, triangles, quadrangles, ...
 *
 * When constant, he can be view as a constant polygon, and thus many existing algorithms can be applied on it.
 * 
 * The main advantage is that points can be accessed with get and set.
 * Therefor, compile-time expression can be used with such objects.
 *
 * His number of points is accessible with the dimension type traits.
 *
 * \remark As it is closed, only unique points have to be specified. Therefore, the begin/end point is no duplicated as for the boost::geometry polygon. In addition, there is no need to go through a ring to access the points.
 *
 */
template < typename Geometry >
class SimpleStaticPolygon
{
    typedef typename std::decay<Geometry>::type polygon_type;
    
    typedef typename point_type<Geometry>::type point_type;

    BOOST_CONCEPT_ASSERT( (concept::Point<point_type>) );
    BOOST_CONCEPT_ASSERT( (boost::RandomAccessRangeConcept<Geometry>) );
    
    enum { pcount = static_num_points<Geometry>::value };   

    /*
    template <
        typename P,
        std::size_t Dimension,
        std::size_t DimensionCount
    >
    struct dimension_checker
    {
        static void apply()
        {
            P* polygon = 0;
            geometry::set<Dimension>(*polygon, geometry::get<Dimension>(*polygon) );
            dimension_checker<P, Dimension+1, DimensionCount>::apply();
        }
    };

    template <
        typename P,
        std::size_t DimensionCount
    >
    struct dimension_checker<P, DimensionCount, DimensionCount>
    {
        static void apply();
    };
    */

    template <
        typename P,
        std::size_t Index,
        std::size_t IndexCount,
        std::size_t Dimension,
        std::size_t DimensionCount
    >
    struct nested_dimension_checker
    {
        static void apply()
        {
            P* polygon = 0;
            geometry::set<Index,Dimension>(*polygon, geometry::get<Index,Dimension>(*polygon) );
            nested_dimension_checker<P, Index, IndexCount, Dimension+1, DimensionCount>::apply();
        }
    };

    template < typename P,
        std::size_t Index,
        std::size_t IndexCount,
        std::size_t DimensionCount
    >
    struct nested_dimension_checker<P, Index, IndexCount, DimensionCount, DimensionCount>
    {
        static void apply()
        {
            nested_dimension_checker<P, Index+1, IndexCount, 0, DimensionCount>::apply();
        }
    };

    template < typename P,
        std::size_t IndexCount,
        std::size_t Dimension,
        std::size_t DimensionCount
    >
    struct nested_dimension_checker<P, IndexCount, IndexCount, Dimension, DimensionCount>
    {
        static void apply();
    };

public:

    BOOST_CONCEPT_USAGE(SimpleStaticPolygon)
    {
        // dimension_checker<Geometry, 0, pcount>::apply();
        nested_dimension_checker < Geometry,
            0, pcount,
            0, dimension<point_type>::value
        >::apply();
    }
};

//! Const version
template < typename Geometry >
class ConstSimpleStaticPolygon
{
    typedef typename std::decay<Geometry>::type polygon_type;
    
    typedef typename point_type<Geometry>::type point_type;

    BOOST_CONCEPT_ASSERT( (concept::ConstPoint<point_type>) );
    BOOST_CONCEPT_ASSERT( (boost::RandomAccessRangeConcept<Geometry>) );
    
    enum { pcount = static_num_points<Geometry>::value };   

    template <
        typename P,
        std::size_t Index,
        std::size_t IndexCount,
        std::size_t Dimension,
        std::size_t DimensionCount
    >
    struct nested_dimension_checker
    {
        static void apply()
        {
            P const* polygon = 0;
            typename coordinate_type<point_type>::type value = geometry::get<Index,Dimension>(*polygon);
            boost::ignore_unused_variable_warning(value);
            
            nested_dimension_checker<P, Index, IndexCount, Dimension+1, DimensionCount>::apply();
        }
    };

    template < typename P,
        std::size_t Index,
        std::size_t IndexCount,
        std::size_t DimensionCount
    >
    struct nested_dimension_checker<P, Index, IndexCount, DimensionCount, DimensionCount>
    {
        static void apply()
        {
            nested_dimension_checker<P, Index+1, IndexCount, 0, DimensionCount>::apply();
        }
    };

    template < typename P,
        std::size_t IndexCount,
        std::size_t Dimension,
        std::size_t DimensionCount
    >
    struct nested_dimension_checker<P, IndexCount, IndexCount, Dimension, DimensionCount>
    {
        static void apply();
    };

public:

    BOOST_CONCEPT_USAGE(ConstSimpleStaticPolygon)
    {
        // dimension_checker<Geometry, 0, pcount>::apply();
        nested_dimension_checker < Geometry,
            0, pcount,
            0, dimension<point_type>::value
        >::apply();
    }
};


}}} // namespace boost::geometry::concept

namespace floe { namespace geometry { namespace concept
{
    using boost::geometry::concept::SimpleStaticPolygon;
    using boost::geometry::concept::ConstSimpleStaticPolygon;

}}} // namespace floe::geometry::concept

#endif // FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_SIMPLE_STATIC_POLYGON_CONCEPT_HPP

