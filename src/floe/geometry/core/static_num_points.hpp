/*!
 * \file floe/geometry/core/static_num_points.hpp
 * \brief Getter for number of points of a static polygon.
 * \author Roland Denis
 */

#ifndef FLOE_GEOMETRY_CORE_STATIC_NUM_POINTS_HPP
#define FLOE_GEOMETRY_CORE_STATIC_NUM_POINTS_HPP

#include <boost/mpl/assert.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/static_assert.hpp>

#include <boost/geometry/util/bare_type.hpp>

namespace boost { namespace geometry
{

namespace traits
{

/*! Traits class indicating the number of points of a static geometry
 */
template <
    typename Geometry,
    typename Enable = void
>
struct static_num_points
{
    BOOST_MPL_ASSERT_MSG
    (
        false, NOT_IMPLEMENTED_FOR_THIS_GEOMETRY_TYPE, (types<Geometry>)
    );
};

} // namespace traits

namespace core_dispatch
{

template <typename T, typename G>
struct static_num_points 
    : static_num_points< simple_static_polygon_tag, G >
{};

template <typename P>
struct static_num_points < simple_static_polygon_tag, P>
    : traits::static_num_points<typename geometry::util::bare_type<P>::type>
{};


} // namespace core_dispatch

template <typename Geometry>
struct static_num_points
    : core_dispatch::static_num_points
        <
            typename tag<Geometry>::type,
            typename geometry::util::bare_type<Geometry>::type
        >
{};

/*! static_num_points, enables compile-time checking if points number are as expected
 */
template <typename Geometry, std::size_t NumPoints>
inline void assert_static_num_points()
{
    BOOST_STATIC_ASSERT
    ((
        boost::mpl::equal_to
        <
            boost::mpl::int_<geometry::static_num_points<Geometry>::value>,
            boost::mpl::int_<NumPoints>
        >::type::value
    ));
}


}} // namespace boost::geometry


namespace floe { namespace geometry
{
    using boost::geometry::static_num_points;
    using boost::geometry::assert_static_num_points;
}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_CORE_STATIC_NUM_POINTS_HPP

