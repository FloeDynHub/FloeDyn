/*!
 * \file floe/geometry/core/point_order.hpp
 * \brief Specify point order (clockwise or counter-clockwise) for circles.
 * \author Roland Denis
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/geometry/reference/core/point_order.html">boost/geometry/core/point_order.hpp</a>
 */

#ifndef FLOE_GEOMETRY_CORE_POINT_ORDER_HPP
#define FLOE_GEOMETRY_CORE_POINT_ORDER_HPP

#include <boost/geometry/core/point_order.hpp>

namespace boost { namespace geometry
{

namespace core_dispatch
{

template <typename TCircle>
struct point_order<circle_tag, TCircle>
    : public detail::point_order::clockwise {};

} // namespace core_dispatch

}} // namespace boost::geometry
#endif // FLOE_GEOMETRY_CORE_POINT_ORDER_HPP

