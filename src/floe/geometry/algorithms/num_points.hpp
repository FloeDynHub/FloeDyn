/*!
 * \file floe/geometry/algorithms/num_points.hpp
 * \brief Return the number of points of a circle, and import other version from boost.
 * \author Roland Denis
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/geometry/reference/algorithms/num_points.html">boost/geometry/algorithms/num_points.hpp</a>
 */

#ifndef FLOE_GEOMETRY_ALGORITHMS_NUM_POINTS_HPP
#define FLOE_GEOMETRY_ALGORITHMS_NUM_POINTS_HPP

#include <boost/geometry/algorithms/num_points.hpp>
#include "floe/geometry/core/tags.hpp"

namespace boost { namespace geometry
{

namespace dispatch
{

template <typename Geometry, bool AddForOpen>
struct num_points<Geometry, AddForOpen, circle_tag>
    : detail::counting::other_count<1>
{};

} // namespace dispatch


}} // namespace boost::geometry


namespace floe { namespace geometry
{

using boost::geometry::num_points;

}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_ALGORITHMS_NUM_POINTS_HPP

