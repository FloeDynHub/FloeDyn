#ifndef FLOE_GEOMETRY_IO_DSV_WRITE_HPP
#define FLOE_GEOMETRY_IO_DSV_WRITE_HPP

#include <ostream>

#include <boost/geometry/io/dsv/write.hpp>

#include "floe/geometry/views/center_view.hpp"
#include "floe/geometry/core/radius.hpp"
#include "floe/geometry/core/tags.hpp"


namespace boost { namespace geometry
{

namespace detail { namespace dsv
{

//! Stream circle as \ref DSV
template < typename TCircle >
struct dsv_circle
{
    template < typename Char, typename Traits >
    static inline
    void apply( 
        std::basic_ostream<Char, Traits>& os,
        TCircle const& circle,
        dsv_settings const& settings
    )
    {
        os << settings.list_open;
        dsv_point< center_view<const TCircle> >::apply( os, center_view<const TCircle>(circle), settings );
        os << settings.list_separator;
        os << get_radius<TCircle>(circle);
        os << settings.list_close;
    }

};

template < typename TSimpleStaticPolygon >
struct dsv_simple_static_polygon
{
    template < typename Char, typename Traits >
    static inline void apply(
        std::basic_ostream<Char, Traits>& os,
        TSimpleStaticPolygon const& polygon,
        dsv_settings const& settings
    )
    {
        dsv_range<TSimpleStaticPolygon>::apply(os, polygon, settings);
    }
};

}} // namespace detail::dsv

namespace dispatch {

template <typename TCircle>
struct dsv<circle_tag, TCircle>
    : detail::dsv::dsv_circle<TCircle>
{};

template <typename TSimpleStaticPolygon>
struct dsv<simple_static_polygon_tag, TSimpleStaticPolygon>
    : detail::dsv::dsv_simple_static_polygon<TSimpleStaticPolygon>
{};

template <typename TTriangle>
struct dsv<triangle_tag, TTriangle>
    : detail::dsv::dsv_simple_static_polygon<TTriangle>
{};

} // namespace dispatch

}} // namespace boost::geometry

namespace floe { namespace geometry
{
    using boost::geometry::dsv;

}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_IO_DSV_WRITE_HPP

