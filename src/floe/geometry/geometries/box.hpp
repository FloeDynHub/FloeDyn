#ifndef FLOE_GEOMETRY_GEOMETRIES_BOX_HPP
#define FLOE_GEOMETRY_GEOMETRIES_BOX_HPP

#include <boost/geometry/geometries/box.hpp>

namespace floe { namespace geometry {

/*! Axis Aligned Box type
 *
 * Use the default box model from boost::geometry
 *
 * \tparam TPoint   Point type
 */

template < typename TPoint >
using Box = boost::geometry::model::box<TPoint>;

}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_GEOMETRIES_BOX_HPP
