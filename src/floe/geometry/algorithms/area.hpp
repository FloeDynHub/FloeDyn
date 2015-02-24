/*!
 * \file floe/geometry/algorithms/area.hpp
 * \brief Area calculation for circles.
 * \author Roland Denis
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/geometry/reference/algorithms/area/area_2_with_strategy.html">boost/geometry/algorithms/area.hpp</a>
 *
 * \namespace boost::geometry::detail::area
 * \brief Implementation details for area calculation.
 */

#ifndef FLOE_GEOMETRY_ALGORITHMS_AREA_HPP
#define FLOE_GEOMETRY_ALGORITHMS_AREA_HPP

#include <boost/geometry/algorithms/area.hpp>
#include <boost/math/constants/constants.hpp>

namespace boost { namespace geometry 
{

namespace detail { namespace area
{

//! Area of a circle
struct circle_area
{
    template <typename Circle, typename Strategy>
    static inline typename coordinate_type<Circle>::type
    apply(Circle const& circle, Strategy const&)
    {
        return boost::math::constants::pi<typename coordinate_type<Circle>::type>() * get_radius(circle) * get_radius(circle);
    }
};

}} // namespace detail::area


namespace dispatch
{

template <typename Geometry>
struct area<Geometry, circle_tag> : detail::area::circle_area
{};

} // namespace dispatch

}} // namespace boost::geometry


namespace floe { namespace geometry
{
    //! Make it visible in our namespace
    using boost::geometry::area;

}} // namespace floe::geometry


#endif // FLOE_GEOMETRY_ALGORITHMS_AREA_HPP
