/*!
 * \file floe/geometry/algorithms/centroid.hpp
 * \brief Centroid for circles and import of boost versions in floe namespace.
 * \author Roland Denis
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/geometry/reference/algorithms/centroid/centroid_3_with_strategy.html">boost/geometry/algorithms/centroid.hpp</a>
 *
 * \namespace boost::geometry::detail::centroid
 * \brief Implementation details for centroid calculation.
 */

#ifndef FLOE_GEOMETRY_ALGORITHM_CENTROID_HPP
#define FLOE_GEOMETRY_ALGORITHM_CENTROID_HPP

#include <boost/geometry/algorithms/centroid.hpp>

#include "floe/geometry/views/center_view.hpp"

namespace boost { namespace geometry {

namespace detail { namespace centroid {

struct centroid_circle
{
    template <
        typename TCircle,
        typename TPointCentroid,
        typename TStrategy
    >
    static inline
    void apply( TCircle const& circle, TPointCentroid & centroid, TStrategy const& )
    {
        geometry::convert( center_view<const TCircle>(circle), centroid );
    }
};

}} // namespace detail::centroid

namespace dispatch
{

template < typename Geometry >
struct centroid< Geometry, circle_tag >
    : detail::centroid::centroid_circle
{};

} // namespace dispatch



}} // namespace boost::geometry

namespace floe { namespace geometry
{
    using boost::geometry::centroid;
    using boost::geometry::return_centroid;

}} // namespace floe::geometry


#endif // FLOE_GEOMETRY_ALGORITHM_CENTROID_HPP
