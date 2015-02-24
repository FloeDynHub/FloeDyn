#ifndef FLOE_GEOMETRY_GEOMETRIES_MULTI_POINT_HPP
#define FLOE_GEOMETRY_GEOMETRIES_MULTI_POINT_HPP

#include <memory>
#include <vector>

#include <boost/geometry/multi/geometries/multi_point.hpp>

namespace floe { namespace geometry
{

template <
    typename TPoint,
    template <typename, typename> class TContainer = std::vector,
    template <typename> class TAllocator = std::allocator
>
using MultiPoint = boost::geometry::model::multi_point<TPoint, TContainer, TAllocator>;

}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_GEOMETRIES_MULTI_POINT_HPP

