/*!
 * \file floe/geometry/algorithms/buffer.hpp
 * \brief Buffer for circles and importation of boost versions in floe namespace.
 * \author Roland Denis
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/geometry/reference/algorithms/buffer/buffer_7_with_strategies.html">boost/geometry/algorithms/buffer.hpp</a>
 *
 * \namespace boost::geometry::detail::buffer
 * \brief Implementation details for buffer calculation.
 */

#ifndef FLOE_GEOMETRY_ALGORITHMS_BUFFER_HPP
#define FLOE_GEOMETRY_ALGORITHMS_BUFFER_HPP

#include <boost/geometry/algorithms/buffer.hpp>

#include "floe/geometry/core/tags.hpp"
#include "floe/geometry/core/radius.hpp"

namespace boost { namespace geometry
{

namespace detail { namespace buffer
{

//! Extends a circle by a specified distance
template <
    typename TCircleIn,
    typename TCircleOut,
    typename T
>
inline
void buffer_circle( TCircleIn const& circle_in, T const& distance, TCircleOut & circle_out )
{
    set<0>(circle_out, get<0>(circle_in));
    set<1>(circle_out, get<1>(circle_in));
    set_radius(circle_out, get_radius(circle_in) + distance);
}


}} // namespace detail::buffer

namespace dispatch {

template <
    typename TCircleIn,
    typename TCircleOut
>
struct buffer< TCircleIn, TCircleOut, circle_tag, circle_tag >
{
    template < typename Distance >
    static inline
    void apply( TCircleIn const& circle_in, Distance const& distance, Distance const&, TCircleOut & circle_out )
    {
        detail::buffer::buffer_circle( circle_in, distance, circle_out );
    }

};

} // namespace dispatch

}} // namespace boost::geometry

namespace floe { namespace geometry
{
    using boost::geometry::buffer;
    using boost::geometry::return_buffer;

}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_ALGORITHMS_BUFFER_HPP
