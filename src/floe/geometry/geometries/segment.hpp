#ifndef FLOE_GEOMETRY_GEOMETRIES_SEGMENT_HPP
#define FLOE_GEOMETRY_GEOMETRIES_SEGMENT_HPP

#include <boost/geometry/geometries/segment.hpp>

namespace floe { namespace geometry
{

    template < typename TPoint >
    using Segment = boost::geometry::model::segment<TPoint>;

}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_GEOMETRIES_SEGMENT_HPP
