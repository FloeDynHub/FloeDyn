/*!
 * \file floe/floes/static_floe.hpp
 * \brief Definition and manipulation of static floe.
 * \author Roland denis
 */

#ifndef FLOE_FLOES_STATIC_FLOE_HPP
#define FLOE_FLOES_STATIC_FLOE_HPP

#include <cmath>

#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/point.hpp"           // Default point type
#include <boost/geometry/geometries/polygon.hpp>        // Default geometry type
#include "floe/geometry/geometries/triangle_mesh.hpp"   // Default mesh type
#include "floe/geometry/frame/theta_frame.hpp"                // Default frame

#include "floe/integration/integrate.hpp"
#include "floe/integration/gauss_legendre.hpp"

namespace floe { namespace floes
{

/*! Static Floe
 *
 * It represent a floe without knowledge of speed, acceleration, etc.
 * The geometry and mesh is defined relatively to its own frame.
 *
 * The origin of the geometry & mesh must be the center of mass (to be cleaned)
 *
 * \tparam T            Fondamental type.
 * \tparam TPoint       Point type.
 * \tparam TGeometry    Geometry type of the floe.
 * \tparam TMesh        Mesh type of the floe.
 * \tparam TFrame       Frame type of the floe.
 * \tparam TDensity     Density type.
 *
 */
template <
    typename T,
    typename TPoint     = geometry::Point<T>,
    typename TGeometry  = boost::geometry::model::polygon<TPoint, false, false>,
    typename TMesh      = geometry::TriangleMesh<TPoint>,
    typename TFrame     = geometry::frame::ThetaFrame<TPoint>,
    typename TDensity   = T
>
class StaticFloe
{

public:

    // Type traits
    typedef T           value_type;
    typedef TPoint      point_type;
    typedef TGeometry   geometry_type;
    typedef TMesh       mesh_type;
    typedef TFrame      frame_type;
    typedef TDensity    density_type;

    using Uptr_geometry_type = std::unique_ptr<geometry_type>;

    //! Default constructor.
    StaticFloe() : m_frame{0,0,0}, m_geometry{nullptr}, m_mesh{nullptr}, m_density{917}, m_mu_static{0.7}, m_area{-1}, m_moment_cst{-1} {}

    //! Deleted copy-constructor.
    StaticFloe( StaticFloe const& ) = delete;

    //! Deleted copy operator.
    StaticFloe& operator= (StaticFloe const) = delete;


    //! Frame accessors
    inline  frame_type const&       get_frame()                             const   { return m_frame; }
    inline  void                    set_frame( frame_type const& frame )            { m_frame = frame; }
    inline  frame_type const&       frame()                                 const   { return m_frame; }
    inline  frame_type &            frame()                                         { return m_frame; }

    //! Geometry accessors
    inline  void                    attach_geometry_ptr( Uptr_geometry_type geometry )  { m_geometry = std::move(geometry); }
    inline  geometry_type const&    geometry()                              const   { return *m_geometry; }
    inline  geometry_type &         geometry()                                      { return *m_geometry; }
    inline  bool                    has_geometry()                          const   { return m_geometry != nullptr; }
    inline  geometry_type const&    get_geometry()                          const   { return *m_geometry; }
    inline  void                    set_geometry( geometry_type const& geometry )
    { 
        if (! has_geometry() ) m_geometry = new geometry_type();
        *m_geometry = geometry;
    }

    //! Mesh accessors
    inline  void                    attach_mesh_ptr( mesh_type* mesh)          { m_mesh = mesh; }
    inline  mesh_type const&        mesh()                      const   { return *m_mesh; }
    inline  mesh_type &             mesh()                              { return *m_mesh; }
    inline  bool                    has_mesh()                  const   { return m_mesh != nullptr; }
    inline  mesh_type const&        get_mesh()                  const   { return *m_mesh; }
    inline  void                    set_mesh( mesh_type& mesh ) 
    { 
        // if (! has_mesh() ) m_mesh = new mesh_type();
        m_mesh = &mesh;
    }

    //! Mu accessors
    inline value_type const&    mu_static() const   { return m_mu_static; }
    inline value_type &         mu_static()         { return m_mu_static; }

    //! Area
    inline value_type area() const { return calc_area(); }

    //! Mass
    inline density_type mass() const { return m_density * area(); }

    //! Momentum constant
    inline value_type moment_cst() const { return calc_moment_cst(); }

    //! Density accessors
    inline value_type   get_density()   const   { return m_density; }
    inline void         set_density( value_type const& density ) { m_density = density; m_moment_cst = -1; }
    inline value_type const&    density()   const { return m_density; }

private:

    frame_type m_frame;         //!< Frame
    Uptr_geometry_type m_geometry;  //!< Geometry (border)
    mesh_type* m_mesh;          //!< Mesh
    density_type m_density;     //!< Density
    value_type m_mu_static;     //!< Static friction coefficient

    mutable value_type m_area;  //!< Area (cached)
    // mutable density_type m_mass; //!< Mass (cached)
    mutable value_type m_moment_cst; //!< Momentum constant (cached)

    //! Calculate area, if not already done.
    inline
    value_type calc_area() const
    {
        return (m_area >= 0) ? m_area : ( has_geometry() ? ( m_area = geometry::area(*m_geometry) ) : m_area );
    }

    //! Calculate momentum constant, if not already done.
    inline
    value_type calc_moment_cst() const
    {
        if ( m_moment_cst >= 0 || ! has_mesh() ) {
            return m_moment_cst;
        }
        else 
        {
            //! \warning To compare results with Matlab code, choose Gauss-Legendre quadrature with 1 point instead of 3 points here.
            const auto strategy = integration::RefGaussLegendre<value_type, 2, 1>(); // Same as the Matlab code.
            //const auto strategy = integration::RefGaussLegendre<value_type, 2, 2>(); // Greater precision.
            const density_type& density = m_density;

            m_moment_cst = integration::integrate(
                    [&density] ( value_type x, value_type y )
                    {
                        return density * ( x*x + y*y );
                    },
                    *m_mesh,
                    strategy
            );
            return m_moment_cst;
        }
    }

};

}} // namespace floe::floes


#endif // FLOE_FLOES_STATIC_FLOE_HPP

