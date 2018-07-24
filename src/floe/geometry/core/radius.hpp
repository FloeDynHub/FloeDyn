/*!
 * \file floe/geometry/core/radius.hpp
 * \brief Radius accessors for circles.
 * \author Roland Denis
 */

#ifndef FLOE_GEOMETRY_CORE_RADIUS_HPP
#define FLOE_GEOMETRY_CORE_RADIUS_HPP

#include <type_traits>
#include <boost/mpl/assert.hpp>

#include <boost/geometry/util/bare_type.hpp>
#include <boost/geometry/core/radius.hpp>

#include "floe/geometry/core/tag.hpp"
#include "floe/geometry/core/tags.hpp"

namespace boost { namespace geometry
{

namespace core_dispatch
{

template < typename Circle >
struct radius_type < circle_tag, Circle >
{
    typedef typename traits::radius_type<Circle>::type type;
};

template < typename Circle, std::size_t Dimension >
struct radius_access < circle_tag, Circle, Dimension, boost::false_type >
    : detail::radius_access<circle_tag, Circle, Dimension>
{
    BOOST_STATIC_ASSERT(Dimension < 2);
};

} // namespace core_dispatch

}} // namespace boost::geometry

namespace floe { namespace geometry
{
    template <typename Geometry>
    inline
    typename boost::geometry::radius_type<Geometry>::type
    get_radius(Geometry const& geometry)
      {
        return boost::geometry::get_radius<0>(geometry);
      }

    template <typename Geometry>
    inline
    void set_radius(Geometry& geometry,
                    typename boost::geometry::radius_type<Geometry>::type const& radius)
      {
        boost::geometry::set_radius<0>(geometry, radius);
      }
}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_CORE_RADIUS_HPP

