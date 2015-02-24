/*!
 * \file floe/geometry/core/access.hpp
 * \brief Accessors for new geometries. Import also boost versions.
 * \author Roland Denis
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/geometry/reference/access.html">boost/geometry/core/access.hpp</a>
 */

#ifndef FLOE_GEOMETRY_CORE_ACCESS_HPP
#define FLOE_GEOMETRY_CORE_ACCESS_HPP

#include <boost/geometry/core/access.hpp>
#include <boost/type_traits/is_pointer.hpp>

#include "floe/geometry/core/tags.hpp"

namespace boost { namespace geometry
{

namespace detail
{

template <
    typename TGeometry,
    typename TCoordinate,
    std::size_t Dimension
>
struct access_non_pointer
{
    static inline TCoordinate get( TGeometry const& geometry )
    {
        return traits::access<TGeometry, Dimension>::get( geometry );
    }

    static inline void set( TGeometry & geometry, TCoordinate const& value )
    {
        traits::access<TGeometry, Dimension>::set( geometry, value );
    }
};


template <
    typename TGeometry,
    typename TCoordinate,
    std::size_t Dimension
>
struct access_pointer
{
    static inline TCoordinate get( TGeometry const* geometry )
    {
        return traits::access<typename boost::remove_pointer<TGeometry>::type, Dimension>::get( *geometry );
    }

    static inline void set( TGeometry* geometry, TCoordinate const& value )
    {
        traits::access<typename boost::remove_pointer<TGeometry>::type, Dimension>::set( *geometry, value );
    }
};


} // namespace detail


// DISPATCH
namespace core_dispatch
{


template <
    typename TCircle,
    typename TCoordinate,
    std::size_t Dimension
>
struct access < circle_tag, TCircle, TCoordinate, Dimension, boost::false_type >
    : detail::access_non_pointer<TCircle, TCoordinate, Dimension>
{};


template <
    typename TCircle,
    typename TCoordinate,
    std::size_t Dimension
>
struct access < circle_tag, TCircle, TCoordinate, Dimension, boost::true_type >
    : detail::access_pointer<TCircle, TCoordinate, Dimension>
{};


template <
    typename TStaticRing,
    typename TCoordinate,
    std::size_t Dimension
>
struct access < static_ring_tag, TStaticRing, TCoordinate, Dimension, boost::false_type >
    : detail::access_non_pointer<TStaticRing, TCoordinate, Dimension>
{};


template <
    typename TStaticRing,
    typename TCoordinate,
    std::size_t Dimension
>
struct access < static_ring_tag, TStaticRing, TCoordinate, Dimension, boost::true_type >
    : detail::access_pointer<TStaticRing, TCoordinate, Dimension>
{};

// That was a bad idea if we want to use it with C array
// => how can we return real[] ??
// Use instead the range concept !!!
/*
template <
    typename TSimpleStaticPolygon,
    typename TCoordinate,
    std::size_t Dimension
>
struct access < simple_static_polygon_tag, TSimpleStaticPolygon, TCoordinate, Dimension, boost::false_type >
    : detail::access_non_pointer<TSimpleStaticPolygon, TCoordinate, Dimension>
{};


template <
    typename TSimpleStaticPolygon,
    typename TCoordinate,
    std::size_t Dimension
>
struct access < simple_static_polygon_tag, TSimpleStaticPolygon, TCoordinate, Dimension, boost::true_type >
    : detail::access_pointer<TSimpleStaticPolygon, TCoordinate, Dimension>
{};
*/

template <
    typename TSimpleStaticPolygon,
    typename TCoordinate,
    std::size_t Index,
    std::size_t Dimension
>
struct indexed_access< simple_static_polygon_tag, TSimpleStaticPolygon, TCoordinate, Index, Dimension, boost::false_type >
    : detail::indexed_access_non_pointer<TSimpleStaticPolygon, TCoordinate, Index, Dimension>
{};


template <
    typename TSimpleStaticPolygon,
    typename TCoordinate,
    std::size_t Index,
    std::size_t Dimension
>
struct indexed_access< simple_static_polygon_tag, TSimpleStaticPolygon, TCoordinate, Index, Dimension, boost::true_type >
    : detail::indexed_access_pointer<TSimpleStaticPolygon, TCoordinate, Index, Dimension>
{};

template <
    typename TSimpleStaticPolygon,
    typename TCoordinate,
    std::size_t Index,
    std::size_t Dimension
>
struct indexed_access< triangle_tag, TSimpleStaticPolygon, TCoordinate, Index, Dimension, boost::false_type >
    : detail::indexed_access_non_pointer<TSimpleStaticPolygon, TCoordinate, Index, Dimension>
{};


template <
    typename TSimpleStaticPolygon,
    typename TCoordinate,
    std::size_t Index,
    std::size_t Dimension
>
struct indexed_access< triangle_tag, TSimpleStaticPolygon, TCoordinate, Index, Dimension, boost::true_type >
    : detail::indexed_access_pointer<TSimpleStaticPolygon, TCoordinate, Index, Dimension>
{};

} // namespace core_dispatch

}} // namespace boost::geometry

namespace floe { namespace geometry
{
    using boost::geometry::get;
    using boost::geometry::set;
}} // namespace floe::geometry


#endif // FLOE_GEOMETRY_CORE_ACCESS_HPP
