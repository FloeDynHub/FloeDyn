/*!
 * \file floe/geometry/algorithms/disjoint.hpp
 * \brief Test if two circles are disjoint, and import other versions from boost.
 * \author Roland Denis
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/geometry/reference/algorithms/disjoint.html">boost/geometry/algorithms/disjoint.hpp</a>
 *
 * \namespace boost::geometry::detail::disjoint
 * \brief Implementation details for disjoint algorithm.
 */

#ifndef FLOE_GEOMETRY_ALGORITHMS_DISJOINT_HPP
#define FLOE_GEOMETRY_ALGORITHMS_DISJOINT_HPP

#include <boost/geometry/algorithms/disjoint.hpp>
#include "floe/geometry/core/tags.hpp"
#include "floe/geometry/views/center_view.hpp"
#include "floe/geometry/core/radius.hpp"

namespace boost { namespace geometry 
{

namespace detail { namespace disjoint
{

template
<
    typename Circle1, typename Circle2,
    std::size_t DimensionCount
>
struct circle_circle
{
    static inline bool apply( Circle1 const& c1, Circle2 const& c2 )
    {
        return 
                boost::geometry::disjoint( center_view<Circle1>(c1), center_view<Circle2>(c2) )
            ||  get_radius(c1) != get_radius(c2);
    }
};

}} // namespace detail::disjoint

namespace dispatch
{

template <typename Circle1, typename Circle2, std::size_t DimensionCount>
struct disjoint<Circle1, Circle2, DimensionCount, circle_tag, circle_tag, false>
    : detail::disjoint::circle_circle<Circle1, Circle2, DimensionCount>
{};

} // namespace dispatch

}} // namespace boost::geometry

#endif // FLOE_GEOMETRY_ALGORITHMS_DISJOINT_HPP

