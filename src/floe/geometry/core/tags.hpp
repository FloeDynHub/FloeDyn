/*!
 * \file floe/geometry/core/tags.hpp
 * \brief Defines tag for new geometries and import those for boost geometries.
 * \author Roland Denis
 * \see boost/geometry/core/tags.hpp
 */

#ifndef FLOE_GEOMETRY_CORE_TAGS_HPP
#define FLOE_GEOMETRY_CORE_TAGS_HPP

#include <boost/geometry/core/tags.hpp>

namespace boost { namespace geometry {

//! Circular tag
struct circular_tag : areal_tag {};

//! Circle tag
struct circle_tag : single_tag, circular_tag {};

//! Multi circle
struct multi_circle_tag : multi_tag, circular_tag {};

//! Single geometry type of multi_circle_tag
template <>
struct single_tag_of<multi_circle_tag>
{
    typedef circle_tag type;
};

//! Static polygon
struct static_polygon_tag : single_tag, polygonal_tag {};

//! Static ring
struct static_ring_tag : single_tag, polygonal_tag {};

//! Multi static polygon
struct multi_static_polygon_tag : multi_tag, polygonal_tag {};

//! Single geometry type of multi_static_polygon_tag
template <>
struct single_tag_of<multi_static_polygon_tag>
{
    typedef static_polygon_tag type;
};

//! Simple static polygon
struct simple_static_polygon_tag : polygon_tag {};

//! Multi simple static polygon
struct multi_simple_static_polygon_tag : multi_tag, polygonal_tag {};

//! Single geometry type of multi_simple_static_polygon_tag
template <>
struct single_tag_of<multi_simple_static_polygon_tag>
{
    typedef simple_static_polygon_tag type;
};

//! Triangle
struct triangle_tag : simple_static_polygon_tag {};

//! Multi triangle
struct multi_triangle_tag : multi_tag, triangle_tag {};

//! Single geometry type of multi_triangle_tag
template <>
struct single_tag_of<multi_triangle_tag>
{
    typedef triangle_tag type;
};

//! Mesh
struct mesh_tag : single_tag, areal_tag {}; 

}} // namespace boost::geometry

namespace floe { namespace geometry
{
    using boost::geometry::circle_tag;
    using boost::geometry::multi_circle_tag;
    using boost::geometry::point_tag;
    using boost::geometry::segment_tag;
    using boost::geometry::box_tag;
    using boost::geometry::polygon_tag;
    using boost::geometry::static_polygon_tag;
    using boost::geometry::static_ring_tag;
    using boost::geometry::triangle_tag;
    using boost::geometry::multi_triangle_tag;
    using boost::geometry::simple_static_polygon_tag;
    using boost::geometry::multi_simple_static_polygon_tag;
    using boost::geometry::mesh_tag;

    // to be continued ...

}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_CORE_TAGS_HPP
