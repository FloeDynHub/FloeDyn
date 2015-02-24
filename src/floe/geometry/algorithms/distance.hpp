/*!
 * \file floe/geometry/algorithms/distance.hpp
 * \brief Distance point/circle, segment/circle and circle/circle. Import algo boost versions.
 * \author Roland Denis
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/geometry/reference/algorithms/distance/distance_3_with_strategy.html">boost/geometry/algorithms/distance.hpp</a>
 *
 * \namespace boost::geometry::detail::distance
 * \brief Implementation details for distance calculation.
 */

#ifndef FLOE_GEOMETRY_ALGORITHMS_DISTANCE_HPP
#define FLOE_GEOMETRY_ALGORITHMS_DISTANCE_HPP

#include "floe/geometry/core/reverse_dispatch.hpp"
#include <boost/geometry/algorithms/distance.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace distance
{

using strategy::distance::services::return_type;

template <typename TPoint, typename TCircle, typename TStrategy>
struct point_to_circle
{
    static inline typename return_type<TStrategy, TPoint, typename point_type<TCircle>::type >::type
    apply( TPoint const& point, TCircle const& circle, TStrategy const& strategy )
    {
        return boost::geometry::distance(point, center_view<TCircle>(circle), strategy) - get_radius(circle);
    }
};

template <typename TCircle1, typename TCircle2, typename TStrategy>
struct circle_to_circle
{
    static inline typename return_type<TStrategy, typename point_type<TCircle1>::type, typename point_type<TCircle2>::type >::type
    apply( TCircle1 const& c1, TCircle2 const& c2, TStrategy const& strategy )
    {
        return boost::geometry::distance(center_view<TCircle1>(c1), center_view<TCircle2>(c2), strategy) - get_radius(c1) - get_radius(c2);
    }
};

template <typename TSegment, typename TCircle, typename TStrategy>
struct segment_to_circle
{
    static inline typename return_type<TStrategy, typename point_type<TSegment>::type, typename point_type<TCircle>::type >::type
    apply( TSegment const& segment, TCircle const& circle, TStrategy const& strategy )
    {
        return boost::geometry::distance(segment, center_view<TCircle>(circle), strategy) - get_radius(circle);
    }
};

}} // namespace detail::distance

namespace dispatch
{

using strategy::distance::services::return_type;

// Point-Circle
template < typename TPoint, typename TCircle, typename TStrategy >
struct distance
    <
        TPoint, TCircle, TStrategy,
        point_tag, circle_tag, strategy_tag_distance_point_point,
        false
    >
    : detail::distance::point_to_circle<TPoint, TCircle, TStrategy>
{};

// Circle-Circle
template < typename TCircle1, typename TCircle2, typename TStrategy >
struct distance
    <
        TCircle1, TCircle2, TStrategy,
        circle_tag, circle_tag, strategy_tag_distance_point_point,
        false
    >
    : detail::distance::circle_to_circle<TCircle1, TCircle2, TStrategy>
{};

// Segment-Circle
template < typename TSegment, typename TCircle, typename TStrategy >
struct distance
    <
        TSegment, TCircle, TStrategy,
        segment_tag, circle_tag, strategy_tag_distance_point_point,
        false
    >
    : detail::distance::segment_to_circle<TSegment, TCircle, TStrategy>
{};

} // namespace dispatch


}} // namespace boost::geometry


namespace floe { namespace geometry
{

using boost::geometry::distance;

}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_ALGORITHMS_DISTANCE_HPP
