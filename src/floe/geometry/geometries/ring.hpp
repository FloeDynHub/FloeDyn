#ifndef FLOE_GEOMETRY_GEOMETRIES_RING_HPP
#define FLOE_GEOMETRY_GEOMETRIES_RING_HPP

#include <boost/geometry/geometries/ring.hpp>

namespace floe { namespace geometry {

/*! Ring type
 *
 * Use the default ring model from boost::geometry
 *
 * \tparam TPoint       Point type
 * \tparam ClockWise    true if the ring is clockwise oriented
 * \tparam Closed       true if the ring is closed
 */
template <
    typename TPoint,
    bool ClockWise = true,
    bool Closed = true
>
using Ring = boost::geometry::model::ring< TPoint, ClockWise, Closed >;

}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_GEOMETRIES_RING_HPP
