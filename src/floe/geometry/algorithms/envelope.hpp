/*!
 * \file floe/geometry/algorithms/envelope.hpp
 * \brief Return the envelope (bounding box) of a circle and import other versions from boost.
 * \author Roland Denis
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/geometry/reference/algorithms/envelope.html">boost/geometry/algorithms/envelope.hpp</a>
 *
 * \namespace boost/geometry/detail/envelope
 * \brief Implementation details for envelope calculations.
 */

#ifndef FLOE_GEOMETRY_ALGORITHMS_ENVELOPE_HPP
#define FLOE_GEOMETRY_ALGORITHMS_ENVELOPE_HPP

#include <boost/geometry/algorithms/envelope.hpp>

#include "floe/geometry/core/tags.hpp"
#include "floe/geometry/core/radius.hpp"

namespace boost { namespace geometry
{

namespace detail { namespace envelope
{

struct envelope_circle
{
    template <
        typename TCircle,
        typename TBox
    >
    static inline
    void apply( TCircle const& circle, TBox & mbr )
    {
        assert_dimension<TBox, 2>();

        typename radius_type<TCircle>::type const radius = get_radius(circle);

        set<min_corner, 0>( mbr, get<0>(circle) - radius );
        set<min_corner, 1>( mbr, get<1>(circle) - radius );
        set<max_corner, 0>( mbr, get<0>(circle) + radius );
        set<max_corner, 1>( mbr, get<1>(circle) + radius );
    }
};

}} // namespace detail::envelope

namespace dispatch {

template <typename TCircle>
struct envelope<TCircle, circle_tag>
    : detail::envelope::envelope_circle
{};

} // namespace dispatch

}} // namespace boost::geometry

namespace floe { namespace geometry
{
    using boost::geometry::envelope;
    using boost::geometry::return_envelope;
    
}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_ALGORITHMS_ENVELOPE_HPP
