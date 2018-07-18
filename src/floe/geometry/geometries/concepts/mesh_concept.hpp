#ifndef FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_MESH_CONCEPT_HPP
#define FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_MESH_CONCEPT_HPP

#include <boost/concept_check.hpp>
#include <boost/type_traits.hpp>

//#include <boost/geometry/iterators/point_iterator.hpp>

#include "floe/geometry/core/cells.hpp"
#include "floe/geometry/core/cells_type.hpp"
#include "floe/geometry/core/point_type.hpp"

namespace boost { namespace geometry { namespace concepts {

/*! Checks mesh concept
 *
 * A mesh is a list of contiguous base shapes that cover a polygon.
 * A mesh model must have:
 * - a point_iterator (TODO)
 * - a cell_iterator
 * - a segment_iterator ???
 * - a boundary, or exterior ring (const polygon) (TODO)
 */
template < typename TMesh >
class Mesh
{

    typedef typename boost::remove_const<TMesh>::type mesh_type;

    typedef typename traits::cells_const_type<mesh_type>::type cells_const_type;
    //typedef typename traits::cells_mutable_type<mesh_type>::type cells_mutable_type;

    typedef typename traits::point_type<mesh_type>::type point_type;

    typedef typename cells_type<TMesh>::type cells_type;

    struct checker
    {
        static inline void apply()
        {
            //mesh_type* mesh = 0;
            mesh_type const* cmesh = 0;

            //cells_mutable_type cells  = traits::cells<TMesh>::get(*mesh);
            cells_const_type   ccells = traits::cells<TMesh>::get(*cmesh);

            //boost::ignore_unused_variable_warning(mesh);
            boost::ignore_unused_variable_warning(cmesh);
            //boost::ignore_unused_variable_warning(cells);
            boost::ignore_unused_variable_warning(ccells);
        }
    };

public:

    BOOST_CONCEPT_USAGE(Mesh)
    {
        checker::apply();
    }
};

template < typename TMesh >
class ConstMesh
{

    typedef typename boost::remove_const<TMesh>::type mesh_type;

    typedef typename traits::cells_const_type<mesh_type>::type cells_const_type;

    typedef typename cells_type<TMesh>::type cells_type;

    struct checker
    {
        static inline void apply()
        {
            mesh_type const* cmesh = 0;

            cells_const_type   ccells = traits::cells<mesh_type>::get(*cmesh);

            boost::ignore_unused_variable_warning(cmesh);
            boost::ignore_unused_variable_warning(ccells);
        }
    };

public:

    BOOST_CONCEPT_USAGE(ConstMesh)
    {
        checker::apply();
    }
};

}}} // namespace boost::geometry::concepts

namespace floe { namespace geometry { namespace concepts
{
    using boost::geometry::concepts::Mesh;
    using boost::geometry::concepts::ConstMesh;
}}} // namespace floe::geometry::concepts

#endif // FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_MESH_CONCEPT_HPP

