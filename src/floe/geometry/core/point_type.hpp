/*!
 * \file floe/geometry/core/point_type.hpp
 * \brief Get point type of new geometries. Also import boost version.
 * \author Roland Denis
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/geometry/reference/core/point_type.html">boost/geometry/core/point_type.hpp</a>
 */

#ifndef FLOE_GEOMETRY_CORE_POINT_TYPE_HPP
#define FLOE_GEOMETRY_CORE_POINT_TYPE_HPP

#include <boost/geometry/core/point_type.hpp>
#include "floe/geometry/core/cells_type.hpp"

namespace boost { namespace geometry { namespace core_dispatch
{

template <typename MultiSimpleStaticPolygon>
struct point_type<multi_simple_static_polygon_tag, MultiSimpleStaticPolygon>
{
    typedef typename boost::geometry::point_type< typename boost::range_value<MultiSimpleStaticPolygon>::type >::type type;
};

template <typename Mesh>
struct point_type < mesh_tag, Mesh >
{
    typedef typename boost::geometry::point_type< typename geometry::cells_type<Mesh>::type >::type type;
};

}}} // namespace boost::geometry::core_dispatch

namespace floe { namespace geometry
{

using boost::geometry::point_type;

}} // floe::geometry

#endif // FLOE_GEOMETRY_CORE_POINT_TYPE_HPP

