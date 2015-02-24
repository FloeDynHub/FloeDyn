/*!
 * \file floe/geometry/core/cells_type.hpp
 * \brief Type traits to get cells type of a mesh.
 * \author Roland Denis
 */

#ifndef FLOE_GEOMETRY_CORE_CELLS_TYPE_HPP
#define FLOE_GEOMETRY_CORE_CELLS_TYPE_HPP

#include <boost/mpl/assert.hpp>
#include <boost/mpl/if.hpp>
#include <boost/range/value_type.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/remove_reference.hpp>

#include "floe/geometry/core/tag.hpp"
#include "floe/geometry/core/tags.hpp"

namespace boost { namespace geometry
{

namespace traits
{

/*! Traits class to indicate the type of a mesh's cells list.
 */
template <typename TGeometry>
struct cells_const_type
{
    BOOST_MPL_ASSERT_MSG
        (
            false, NOT_IMPLEMENTED_FOR_THIS_GEOMETRY_TYPE
            , (types<TGeometry>)
        );
};

template <typename TGeometry>
struct cells_mutable_type
{
    BOOST_MPL_ASSERT_MSG
        (
            false, NOT_IMPLEMENTED_FOR_THIS_GEOMETRY_TYPE
            , (types<TGeometry>)
        );
};

} // namespace traits


namespace core_dispatch
{

template <typename GeometryTag, typename TGeometry>
struct cells_return_type
{};

template <typename TMesh>
struct cells_return_type<mesh_tag, TMesh>
{
    typedef typename boost::remove_const<TMesh>::type nc_mesh_type;

    typedef typename mpl::if_
        <
            boost::is_const<TMesh>,
            typename traits::cells_const_type<nc_mesh_type>::type,
            typename traits::cells_mutable_type<nc_mesh_type>::type
        >::type type;
};

template <typename GeometryTag, typename TGeometry>
struct cells_type
{};

template <typename TMesh>
struct cells_type<mesh_tag, TMesh>
{
    typedef typename boost::remove_reference
        <
            typename cells_return_type<mesh_tag, TMesh>::type
        >::type type;
};


} // namespace core_dispatch


/*! A mesh is composed of a list of points and cells.
 * This metafunction retrieves the type of the cells list.
 */
template <typename TGeometry>
struct cells_type
{
    typedef typename core_dispatch::cells_type
        <
            typename tag<TGeometry>::type,
            TGeometry
        >::type type;
};

template <typename TGeometry>
struct cells_return_type
{
    typedef typename core_dispatch::cells_return_type
        <
            typename tag<TGeometry>::type,
            TGeometry
        >::type type;
};

}} // namespace boost::geometry


namespace floe { namespace geometry 
{
    using boost::geometry::cells_type;
    using boost::geometry::cells_return_type;
}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_CORE_CELLS_TYPE_HPP

