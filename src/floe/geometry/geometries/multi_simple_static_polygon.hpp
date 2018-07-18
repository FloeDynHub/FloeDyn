#ifndef FLOE_GEOMETRY_GEOMETRIES_MULTI_SIMPLE_STATIC_POLYGON_HPP
#define FLOE_GEOMETRY_GEOMETRIES_MULTI_SIMPLE_STATIC_POLYGON_HPP

#include <memory>
#include <vector>

#include <boost/concept/assert.hpp>

#include "floe/geometry/core/tags.hpp"
#include "floe/geometry/geometries/concepts/multi_simple_static_polygon_concept.hpp"

namespace floe { namespace geometry 
{


/*! MultiSimpleStaticPolygon, a collection a SimpleStaticPolygon
 */
template <
    typename Polygon,
    template <typename, typename> class Container = std::vector,
    template <typename> class Allocator = std::allocator
>
struct MultiSimpleStaticPolygon
    : public Container<Polygon, Allocator<Polygon> >
{
    BOOST_CONCEPT_ASSERT( (concepts::SimpleStaticPolygon<Polygon>) );    
};


}}; // namespace floe::geometry

namespace boost { namespace geometry { namespace traits
{

//! Tag
template <
    typename Polygon,
    template <typename, typename> class Container,
    template <typename> class Allocator
>
struct tag< floe::geometry::MultiSimpleStaticPolygon<Polygon, Container, Allocator> >
{
    typedef multi_simple_static_polygon_tag type;
};

}}} // namespace boost::geometry::traits

#endif // FLOE_GEOMETRY_GEOMETRIES_MULTI_SIMPLE_STATIC_POLYGON_HPP

