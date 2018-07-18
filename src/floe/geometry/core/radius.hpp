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

#include "floe/geometry/core/tag.hpp"
#include "floe/geometry/core/tags.hpp"

namespace boost { namespace geometry
{

namespace traits
{

template <typename Geometry, typename T>
struct radius_access
{
    BOOST_MPL_ASSERT_MSG
    (
        false, NOT_IMPLEMENTED_FOR_THIS_CIRCLE_TYPE, (types<Geometry>)
    );
};

template <typename Geometry>
struct radius_type
{
    BOOST_MPL_ASSERT_MSG
    (
        false, NOT_IMPLEMENTED_FOR_THIS_CIRCLE_TYPE, (types<Geometry>)
    );
};

} // namespace traits

namespace detail 
{

template < typename Circle, typename T >
struct radius_access_non_pointer
{
    static inline T get(Circle const& c)
    {
        return traits::radius_access<Circle, T>::get(c);
    }

    static inline void set( Circle & c, T const& radius )
    {
        traits::radius_access<Circle, T>::set(c, radius);
    }
};


template < typename Circle, typename T >
struct radius_access_pointer
{
    static inline T get(Circle const* c)
    {
        return traits::radius_access<typename std::remove_pointer<Circle>::type, T>::get(*c);
    }

    static inline void set( Circle* c, T const& radius )
    {
        traits::radius_access<typename std::remove_pointer<Circle>::type, T>::set(*c, radius);
    }
};

} // namespace detail

namespace core_dispatch
{

template < typename Tag, typename Geometry >
struct radius_type
{};

template < typename Tag, typename Geometry, typename T, bool IsPointer >
struct radius_access
{};

template < typename Circle >
struct radius_type < circle_tag, Circle >
{
    typedef typename traits::radius_type<Circle>::type type;
};

template < typename Circle, typename T >
struct radius_access < circle_tag, Circle, T, false >
    : detail::radius_access_non_pointer<Circle, T>
{};

template < typename Circle, typename T >
struct radius_access < circle_tag, Circle, T, true >
    : detail::radius_access_pointer<Circle, T>
{};

/*
template < typename Circle, typename T >
struct radius_access < circle_tag, Circle, T >
{
    static inline T get(Circle const& c)
    {
        return traits::radius_access<Circle, T>::get(c);
    }

    static inline void set( Circle & c, T const& radius )
    {
        traits::radius_access<Circle, T>::set(c, radius);
    }
};
*/

} // namespace core_dispatch

template < typename Geometry >
struct radius_type
{
    typedef 
        typename core_dispatch::radius_type < 
            typename tag<Geometry>::type,
            typename boost::geometry::util::bare_type<Geometry>::type
        >::type type;
};

template < typename Geometry >
inline
typename radius_type<Geometry>::type
    get_radius( Geometry const& geometry )
{
    return core_dispatch::radius_access <
        typename tag<Geometry>::type,
        typename boost::geometry::util::bare_type<Geometry>::type,
        typename radius_type<Geometry>::type,
        std::is_pointer<Geometry>::value
    >::get( geometry );
}

template < typename Geometry >
inline
void  set_radius( 
    Geometry & geometry,
    typename radius_type<Geometry>::type const& radius
)
{
    core_dispatch::radius_access <
        typename tag<Geometry>::type,
        typename boost::geometry::util::bare_type<Geometry>::type,
        typename radius_type<Geometry>::type,
        std::is_pointer<Geometry>::value
    >::set( geometry, radius );
}

}} // namespace boost::geometry

namespace floe { namespace geometry
{
    using boost::geometry::get_radius;
    using boost::geometry::set_radius;
}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_CORE_RADIUS_HPP

