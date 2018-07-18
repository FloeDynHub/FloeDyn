#ifndef FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_STATIC_POLYGON_CONCEPT_HPP
#define FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_STATIC_POLYGON_CONCEPT_HPP

#include <cstddef>
#include <type_traits>

#include <boost/concept_check.hpp>

#include <boost/geometry/core/point_type.hpp>

#include "floe/geometry/core/access.hpp"
#include "floe/geometry/core/exterior_ring.hpp"
#include "floe/geometry/core/ring_type.hpp"

#include "floe/geometry/geometries/concepts/static_ring_concept.hpp"
#include <boost/geometry/geometries/concepts/point_concept.hpp>

namespace boost { namespace geometry { namespace concept {

/*! Checks static polygon concept (const version)
 *
 * A static polygon is a polygon made of a fixed number of points
 * and composed only by an exterior static ring.
 * It can be, for example, triangles, quadrangles, ...
 * 
 * \remark We can imagine adding interior ring but there is no use to have
 * only fixed size interior rings, with same size as the interior one.
 *
 * When constant, he can be view as a constant polygon, and thus many existing algorithms can be applied on it.
 *
 * \remark As for the boost::geometry polygon, the first and last point must be identical if the polygon is closed.
 * Therefore, a triangle is, for example, composed of 4 points.
 */
template < typename Geometry >
class StaticPolygon
{
    typedef typename std::decay<Geometry>::type spolygon_type;
    
    // Exterior ring
    typedef typename traits::ring_const_type<spolygon_type>::type sring_const_type;
    typedef typename traits::ring_mutable_type<spolygon_type>::type sring_mutable_type;
    
    typedef typename point_type<Geometry>::type point_type;
    typedef typename ring_type<Geometry>::type sring_type;

    BOOST_CONCEPT_ASSERT( (concept::Point<point_type>) );
    BOOST_CONCEPT_ASSERT( (concept::StaticRing<sring_type>) );

    struct checker
    {
        static inline void apply()
        {
            spolygon_type* poly = 0;
            spolygon_type const* cpoly = poly;

            sring_type se        = traits::exterior_ring<Geometry>::get(*poly);
            sring_const_type sce = traits::exterior_ring<Geometry>::get(*cpoly);

            boost::ignore_unused_variable_warning(se);
            boost::ignore_unused_variable_warning(sce);
        }
    };
    
public:

    BOOST_CONCEPT_USAGE(StaticPolygon)
    {
        checker::apply();
    }
};


template < typename Geometry >
class ConstStaticPolygon
{
    typedef typename std::decay<Geometry>::type cspolygon_type;
    
    // Exterior ring
    typedef typename traits::ring_const_type<cspolygon_type>::type sring_const_type;
    
    typedef typename point_type<cspolygon_type>::type point_type;
    typedef typename ring_type<cspolygon_type>::type sring_type;

    BOOST_CONCEPT_ASSERT( (concept::ConstPoint<point_type>) );
    BOOST_CONCEPT_ASSERT( (concept::ConstStaticRing<sring_type>) );

    struct checker
    {
        static inline void apply()
        {
            cspolygon_type const* cpoly = 0;

            sring_const_type sce = traits::exterior_ring<cspolygon_type>::get(*cpoly);

            boost::ignore_unused_variable_warning(sce);
        }
    };
    
public:

    BOOST_CONCEPT_USAGE(ConstStaticPolygon)
    {
        checker::apply();
    }
};

}}} // namespace boost::geometry::concepts

#endif // FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_STATIC_POLYGON_CONCEPT_HPP
