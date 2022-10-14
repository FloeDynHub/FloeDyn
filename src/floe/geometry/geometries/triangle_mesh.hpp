#ifndef FLOE_GEOMETRY_GEOMETRIES_TRIANGLE_MESH_HPP
#define FLOE_GEOMETRY_GEOMETRIES_TRIANGLE_MESH_HPP

#include <cstddef>
#include <vector>
#include <array>

#include "floe/geometry/core/tag.hpp"

#include <boost/geometry/geometries/multi_point.hpp>
#include "floe/geometry/geometries/multi_simple_static_polygon.hpp"
#include "floe/geometry/geometries/multi_ssp_point_cloud.hpp"
#include "floe/geometry/geometries/triangle.hpp"

#include "floe/geometry/core/cells_type.hpp"
#include "floe/geometry/core/cells.hpp"
#include "floe/geometry/core/point_type.hpp"

namespace floe { namespace geometry {

/*! Mesh model composed of triangles
 *
 * It is a points list with a connectivity list.
 *
 * \tparam TPoint   Point type
 *
 * \todo new triangles must be tested for clockwise orientation
 * \todo mutable cells list that autorize to modify underlying points
 */
template <
    typename TPoint
>
struct TriangleMesh
{
public:
   
    //! Type traits
    using point_type       = TPoint;
    using multi_point_type = typename boost::geometry::model::multi_point<point_type>;
    using cell_type        = typename floe::geometry::Triangle<point_type, true>;
    using connectivity_type = std::vector< std::array<std::size_t,3> >;
    using multi_cell_type  = typename floe::geometry::MultiSSPPointCloud<cell_type, multi_point_type const*, connectivity_type const*>;

    TriangleMesh() : m_points{}, m_connect{}, m_cells{nullptr, nullptr} {}

    //! Return the points list
    inline multi_point_type&       points()       { return m_points; }
    inline multi_point_type const& points() const { return m_points; }

    /*! Return a iterable list of cells
     *
     * \todo please, find another solution ... this is horrible. Make a mutable cells container ?
     */
    inline multi_cell_type const& cells() const
    {
        m_cells = multi_cell_type{ &m_points, &m_connect };
        return m_cells;
    }

    //! Add a triangle (no orientation check)
    inline
    TriangleMesh<TPoint> &
        add_triangle( std::size_t i, std::size_t j, std::size_t k )
    {
        m_connect.push_back( {{i, j, k}} );
        return *this;
    }

    //! Return the connectivity list
    inline connectivity_type&       connectivity()          { return m_connect; }
    inline connectivity_type const& connectivity() const    { return m_connect; }

    inline std::size_t const& size() const    { return m_cells.size(); }

private:

    multi_point_type m_points;
    connectivity_type m_connect;
    mutable multi_cell_type m_cells;
};

}} // namespace floe::geometry

// Type traits
namespace {
    using floe::geometry::TriangleMesh;
}

namespace boost { namespace geometry { namespace traits
{

template <typename TPoint>
struct tag< TriangleMesh<TPoint> >
{
    typedef mesh_tag type;
};

template <typename TPoint>
struct cells_const_type< TriangleMesh<TPoint> >
{
    typedef 
        typename TriangleMesh<TPoint>::multi_cell_type const& 
        type;
};

template <typename TPoint>
struct cells_mutable_type< TriangleMesh<TPoint> >
{
    typedef
        typename TriangleMesh<TPoint>::multi_cell_type const&
        type;
};

template <typename TPoint>
struct cells< TriangleMesh<TPoint> >
{
    typedef TriangleMesh<TPoint> mesh_type;

    static inline
    typename mesh_type::multi_cell_type const& get ( mesh_type& mesh )
    {
        return mesh.cells();
    }

    static inline
    typename mesh_type::multi_cell_type const&
        get ( mesh_type const& mesh )
    {
        return mesh.cells();
    }
};

template <typename TPoint>
struct point_type< TriangleMesh<TPoint> >
{
    typedef TPoint type;
};

}}} // namespace boost::geometry


#endif //FLOE_GEOMETRY_GEOMETRIES_MESH_HPP

