#ifndef FLOE_GEOMETRY_GEOMETRIES_TRIANGLE_MESH_HPP
#define FLOE_GEOMETRY_GEOMETRIES_TRIANGLE_MESH_HPP

#include <cstddef>
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Sparse>



#include "floe/geometry/core/tag.hpp"

#include <boost/geometry/geometries/multi_point.hpp>
#include "floe/geometry/geometries/multi_simple_static_polygon.hpp"
#include "floe/geometry/geometries/multi_ssp_point_cloud.hpp"
#include "floe/geometry/geometries/triangle.hpp"

#include "floe/geometry/core/cells_type.hpp"
#include "floe/geometry/core/cells.hpp"
#include "floe/geometry/core/point_type.hpp"
#include "floe/integration/gauss_legendre.hpp"

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
    using real_type = typename TPoint::value_type; 
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

    // integration functions 
    inline Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> get_integrationMatrix() const {return m_integrationMatrix;};    
    inline Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> get_jacobians() const {return m_jacobians;};    
    inline real_type get_jacobian(size_t iElem) const 
    {
        if (iElem < m_connect.size())
            return triangle_det_jac(iElem);
        else 
            return 0.0 ;
    };    
    inline Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> get_weightAndJac() const {return m_weightAndJac;};    
    inline Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> get_shapeFunAtGP() const {return m_shapeFunAtGP;};    
    bool prepare()
    {
        floe::integration::RefGaussLegendre<double,2,2> i;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> trianglePoints = i.pointsAndWeights();
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> N;
        
        size_t nIntegrationPoints(trianglePoints.rows());
        m_integrationPointsAndWeights.resize(m_connect.size()*nIntegrationPoints, 3);
        m_integrationMatrix.resize(1, m_connect.size()*3);
        m_jacobians.resize(m_connect.size(), 1);
        N.resize(nIntegrationPoints, 3);
        m_shapeFunAtGP.resize(m_connect.size()*nIntegrationPoints,3);
        m_weightAndJac.resize(1, m_connect.size()*nIntegrationPoints);

        for(size_t iElem = 0 ; iElem < m_connect.size() ; ++iElem)
        {
            N.setZero();
            m_integrationPointsAndWeights.block((iElem)*nIntegrationPoints, 0, nIntegrationPoints, 3) = trianglePoints;
            m_jacobians(iElem, 0) = triangle_det_jac(iElem);
            N = buildN(m_integrationPointsAndWeights.block((iElem)*nIntegrationPoints, 0, nIntegrationPoints, 2));
            m_integrationMatrix.block(0, (iElem)*nIntegrationPoints, 1, 3) = m_jacobians(iElem, 0)*m_integrationPointsAndWeights.block((iElem)*nIntegrationPoints, 2, nIntegrationPoints, 1).transpose()*N;
            // the integration vector for one element contains the jacobian, the weights, the fem shape functions. It can be multiplied by the nodal values of a function to integrate it 
            m_shapeFunAtGP.block((iElem)*nIntegrationPoints, 0, nIntegrationPoints, 3) = N;
            m_weightAndJac.block(0, (iElem)*nIntegrationPoints, 1, nIntegrationPoints) = m_jacobians(iElem, 0)*trianglePoints.block(0,2, 3,1).transpose();
        }
        return true;
    }
    //! Return the mesh sizes
    inline size_t get_n_cells() const { return m_connect.size(); }
    inline size_t get_n_nodes() const { return boost::geometry::num_points(m_points); }

private:

    multi_point_type m_points;
    connectivity_type m_connect;
    mutable multi_cell_type m_cells;

    // il faudrait faire rentrer le template real_type 
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_integrationPointsAndWeights;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_shapeFunAtGP; 
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_weightAndJac; 
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_integrationMatrix;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_jacobians;

    
    real_type triangle_det_jac(size_t iElem) const
    {
        double x0 = m_points[m_connect[iElem][0]][0];
        double x1 = m_points[m_connect[iElem][1]][0];
        double x2 = m_points[m_connect[iElem][2]][0];
        double y0 = m_points[m_connect[iElem][0]][1];
        double y1 = m_points[m_connect[iElem][1]][1];
        double y2 = m_points[m_connect[iElem][2]][1];
        return ((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0));
    }

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> buildN(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> points)
    {
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> N;
        size_t nPoints = points.rows();
        N.resize(nPoints, 3);
        if (points.cols() !=2)
            std::cout << "Unexpected matrix size, line " << __LINE__ << " in file " << __FILE__ << ". I got " << points.cols() << " instead of 2" << std::endl;

        for(size_t iPoint = 0; iPoint < nPoints ; ++iPoint)
        {
            N(iPoint, 0) = 1 - points(iPoint, 0) - points(iPoint, 1);
            N(iPoint, 1) = points(iPoint, 0);
            N(iPoint, 2) = points(iPoint, 1);
        }
        return N;
    }
    // size_t m_n_cells;
    // size_t m_n_nodes;
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

