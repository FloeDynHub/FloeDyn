/*!
 * \file floe/geometry/core/cells.hpp
 * \brief Cells accessors for meshes.
 * \author Roland Denis
 */

#ifndef FLOE_GEOMETRY_CORE_CELLS_HPP
#define FLOE_GEOMETRY_CORE_CELLS_HPP

#include <boost/mpl/assert.hpp>
#include <boost/type_traits/remove_const.hpp>

#include <boost/geometry/util/add_const_if_c.hpp>

#include "floe/geometry/core/cells_type.hpp"
#include "floe/geometry/core/tag.hpp"
#include "floe/geometry/core/tags.hpp"

namespace boost { namespace geometry
{

namespace traits
{

/*! Traits class defining access to cell's list of a mesh
 */
template <typename TMesh>
struct cells
{
    BOOST_MPL_ASSERT_MSG
    (
        false, NOT_IMPLEMENTED_FOR_THIS_MESH_TYPE, (types<TMesh>)
    );
};

} // namespace traits

namespace core_dispatch
{

template <typename Tag, typename TMesh>
struct cells
{
    BOOST_MPL_ASSERT_MSG
    (
        false, NOT_IMPLEMENTED_FOR_THIS_MESH_TYPE, (types<TMesh>)
    );
};

template <typename TMesh>
struct cells<mesh_tag, TMesh>
{
    static
    typename geometry::cells_return_type<TMesh>::type
    apply( typename add_const_if_c
        <
            boost::is_const<TMesh>::type::value,
            TMesh
        >::type& mesh )
    {
        return traits::cells
            <
                typename boost::remove_const<TMesh>::type
            >::get(mesh);
    }
};

} // namespace core_dispatch


/*! Function to get the cell's list of a mesh
 *
 * \tparam TMesh mesh type
 * \param mesh the mesh to get the cell's list from
 * \return a reference to a model of multi_simple_static_polygon
 */
template <typename TMesh>
inline typename cells_return_type<TMesh>::type
    cells( TMesh& mesh )
{
    return core_dispatch::cells
        <
            typename tag<TMesh>::type,
            TMesh
        >::apply(mesh);
}

/*! Function to get the cell's list of a mesh (const version)
 *
 * \tparam TMesh mesh type
 * \param mesh the mesh to get the cell's list from
 * \return a reference to a model of multi_simple_static_polygon
 */
template <typename TMesh>
inline typename cells_return_type<TMesh const>::type
    cells( TMesh const& mesh )
{
    return core_dispatch::cells
        <
            typename tag<TMesh>::type,
            TMesh const
        >::apply(mesh);
}


}} //namespace boost::geometry

namespace floe { namespace geometry
{
    using boost::geometry::cells;
}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_CORE_CELLS_HPP

