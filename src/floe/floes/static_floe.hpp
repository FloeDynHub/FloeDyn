/*!
 * \file floe/floes/static_floe.hpp
 * \brief Definition and manipulation of static floe.
 * \author Roland denis
 */

#ifndef FLOE_FLOES_STATIC_FLOE_HPP
#define FLOE_FLOES_STATIC_FLOE_HPP

#include <cmath>
#include <iostream>
#include <vector>

#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/point.hpp"           // Default point type
#include <boost/geometry/geometries/polygon.hpp>        // Default geometry type
#include "floe/geometry/geometries/triangle_mesh.hpp"   // Default mesh type
#include "floe/geometry/frame/theta_frame.hpp"                // Default frame

#include "floe/integration/integrate.hpp"
#include "floe/integration/gauss_legendre.hpp"
#include "floe/generator/mesh_generator.hpp"
#include "floe/geometry/arithmetic/point_operators.hpp"

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
    typedef T           real_type;
    typedef TPoint      point_type;
    typedef TGeometry   geometry_type;
    typedef TMesh       mesh_type;
    typedef TFrame      frame_type;
    typedef TDensity    density_type;

    using Uptr_geometry_type = std::unique_ptr<geometry_type>;

    //! Default constructor.
    StaticFloe() : m_frame{0,0,0}, m_geometry{nullptr}, m_mesh{nullptr}, m_density{917}, m_mu_static{0.7},
                    m_thickness{1}, m_C_w{5 * 1e-3}, m_area{-1}, m_moment_cst{-1}, m_min_crack_energy{0} {}

    StaticFloe(geometry_type new_border) : m_frame{0,0,0}, m_geometry{new_border}, m_mesh{nullptr}, m_density{917}, m_mu_static{0.7},
                    m_thickness{1}, m_C_w{5 * 1e-3}, m_area{-1}, m_moment_cst{-1}, m_min_crack_energy{0} {}


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
    inline  geometry_type &         geometry()                                      { this->reset_area(); return *m_geometry; }
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
    inline real_type const&    mu_static() const   { return m_mu_static; }
    inline real_type &         mu_static()         { return m_mu_static; }
    inline void set_mu_static(real_type mu_static) { m_mu_static = mu_static; }

    //! Thickness accessors
    inline real_type thickness() const { return m_thickness; }
    inline void set_thickness(real_type v){ m_thickness = v; }

    //! Oceanic skin drag accessors
    inline real_type C_w() const { return m_C_w; }
    inline void set_C_w(real_type v){ m_C_w = v; }

    //! Area
    inline real_type area() const { return (m_area >= 0) ? m_area : calc_area(); }

    //! Mass
    inline density_type mass() const { return m_density * m_thickness * area(); }

    //! Momentum constant
    inline real_type moment_cst() const { return calc_moment_cst(); }

    //! Density accessors
    inline real_type   get_density()   const   { return m_density; }
    inline void         set_density( real_type const& density ) { m_density = density; m_moment_cst = -1; }
    inline real_type const&    density()   const { return m_density; }

    inline void update_caracteristic(real_type init_density,real_type init_mu_static, real_type init_thickness, real_type init_C_w);
    // inline void generate_mesh();
    inline point_type get_mass_center() ;
    inline std::vector<TGeometry> fracture_floe();

    real_type min_radius() const
    {
        real_type min_radius = std::numeric_limits<real_type>::max();
        for (size_t i = 0; i < m_geometry->outer().size(); ++i)
        {
            auto const& p = m_geometry->outer()[i];
            for (size_t j = 0; j < m_geometry->outer().size(); ++j) {
                auto const& q = m_geometry->outer()[j];
                if ((i - j) % m_geometry->outer().size() < m_geometry->outer().size() / 3) // TODO : improve this logic
                    continue; // Too close to be a valid radius
                min_radius = std::min(min_radius, norm2(p - q));
            }
        }
        return min_radius;
    }

    //! Area
    inline real_type min_crack_energy() const { return (m_min_crack_energy > 0) ? m_min_crack_energy : calc_min_crack_energy(); }

private:

    frame_type m_frame;         //!< Frame
    Uptr_geometry_type m_geometry;  //!< Geometry (border)
    mesh_type* m_mesh;          //!< Mesh
    density_type m_density;     //!< Density
    real_type m_mu_static;     //!< Static friction coefficient
    real_type m_thickness;     //!< Vertical thickness (constant over floe surface)
    real_type m_C_w;     //!< Oceanic skin drag coeff

    mutable real_type m_area;  //!< Area (cached)
    // mutable density_type m_mass; //!< Mass (cached)
    mutable real_type m_moment_cst; //!< Momentum constant (cached)
    mutable real_type m_min_crack_energy; //!< Minimum crack energy

    //! Calculate area, if not already done.
    inline
    real_type calc_area() const
    {
        return ( has_geometry() ? ( m_area = geometry::area(*m_geometry) ) : m_area );
    }
    inline void reset_area() { m_area = -1; }

    //! Calculate momentum constant, if not already done.
    inline
    real_type calc_moment_cst() const
    {
        if ( m_moment_cst >= 0 || ! has_mesh() ) {
            return m_moment_cst;
        }
        else
        {
            //! \warning To compare results with Matlab code, choose Gauss-Legendre quadrature with 1 point instead of 3 points here.
            const auto strategy = integration::RefGaussLegendre<real_type, 2, 1>(); // Same as the Matlab code.
            //const auto strategy = integration::RefGaussLegendre<real_type, 2, 2>(); // Greater precision.
            const density_type& density = m_density;

            m_moment_cst = integration::integrate(
                    [&density] ( real_type x, real_type y )
                    {
                        return density * ( x*x + y*y );
                    },
                    *m_mesh,
                    strategy
            );
            return m_moment_cst;
        }
    }

    //! Calculate min crack energy
    real_type calc_min_crack_energy() const
    {   
        // Minimum crack energy
        // Source coeff : MODÉLISATION DE LA FRACTURE DE LA GLACE DE MER PAR LA HOULE, Alexandre TLILI, 2022
        real_type ice_crack_coeff = 1.5 * 1000; // Gc entre 1,5 J m−2 et 3,5 J m−2
        m_min_crack_energy = ice_crack_coeff * m_thickness * this->min_radius();
        return m_min_crack_energy;
    }


};


template <typename T,typename TPoint,typename TGeometry,typename TMesh,typename TFrame ,typename TDensity>
std::vector<TGeometry>
StaticFloe<T,TPoint,TGeometry,TMesh,TFrame,TDensity>::fracture_floe()
{
    // Better basic fracture : cutting floe according to crack geometry
    auto& boundary = this->geometry().outer();
    point_type middle_point = (boundary[0] + boundary[boundary.size() - 1]) / 2;
    real_type min_dist = norm2(middle_point);
    point_type crack_start = middle_point;
    // crack_start will be the closest edge midpoint to floe's mass center ({0, 0})
    for (std::size_t i = 0; i < this->geometry().outer().size() - 1; ++i){
        middle_point = (boundary[i] + boundary[i + 1])  / 2;
        if (norm2(middle_point) < min_dist) {
            min_dist = norm2(middle_point);
            crack_start = middle_point;
        }
    }
    std::vector<TGeometry> new_borders;
    // crack is a long and thin rectangle containing crack_start and floe's mass center
    geometry_type crack;
    real_type crack_width = std::sqrt(this->area()) * 0.002;
    point_type crack_ortho = direct_orthogonal(crack_start) / norm2(crack_start);
    point_type crack_delta = crack_ortho * crack_width / 2;
    crack.outer().push_back(crack_start * 1e6 + crack_delta);
    crack.outer().push_back(crack_start * 1e6 - crack_delta);
    crack.outer().push_back(- crack_start * 1e6 - crack_delta);
    crack.outer().push_back(- crack_start * 1e6 + crack_delta);
    boost::geometry::correct(crack);
    // remove crack from floe geometry
    boost::geometry::difference(this->geometry().outer(), crack, new_borders);
    return new_borders;
}


template <typename T,typename TPoint,typename TGeometry,typename TMesh,typename TFrame ,typename TDensity>
void
StaticFloe<T,TPoint,TGeometry,TMesh,TFrame,TDensity>::update_caracteristic(real_type init_density,real_type init_mu_static, real_type init_thickness, real_type init_C_w){

    // generate mesh : récup dans mesh generator
    // this->generate_mesh(); // TODO : might be a good idea !
    this->set_mu_static(init_mu_static);
    this->set_thickness(init_thickness);
    this->set_density( init_density);
    this->set_C_w(init_C_w);
    this->m_area=this->area();
    this->m_moment_cst=this->moment_cst();
}



// template <typename T,typename TPoint,typename TGeometry,typename TMesh,typename TFrame ,typename TDensity>
// void
// StaticFloe<T,TPoint,TGeometry,TMesh,TFrame,TDensity>::generate_mesh()
// {
//         // suppose here border have around 25 vertex
//         TMesh mesh =floe::generator::generate_mesh_for_shape<TGeometry,TMesh>(this->geometry());
//         using integration_strategy = floe::integration::RefGaussLegendre<real_type,2,2>;
//         // Center mesh and shape on floe's center of mass
//         TPoint mass_center=this->get_mass_center();
//         // a faire !!
//         //TGeometry shape_cpy = m_geometry();
//         //Tframe new_frame  {-mass_center, 0};
//         //geometry::transform( *shape_cpy, this->geometry(), geometry::frame::transformer( newframe));
//         //TMesh mesh_cpy = mesh;
//         //geometry::transform( mesh_cpy, mesh, geometry::frame::transformer( typename floe_type::frame_type{-mass_center, 0} ));
//         this->set_mesh(mesh);
// }





template <typename T,typename TPoint,typename TGeometry,typename TMesh,typename TFrame ,typename TDensity>
TPoint
StaticFloe<T,TPoint,TGeometry,TMesh,TFrame,TDensity>::get_mass_center()
{
        using integration_strategy = floe::integration::RefGaussLegendre<real_type,2,2>;
        TPoint mass_center = floe::integration::integrate(
            [] (real_type x, real_type y) { return point_type{x, y}; },
            m_mesh, integration_strategy()
        ) / floe::integration::integrate(
            [] (real_type x, real_type y) { return 1.; },
            m_mesh,integration_strategy()
        );
        return mass_center;
}




}} // namespace floe::floes


#endif // FLOE_FLOES_STATIC_FLOE_HPP

