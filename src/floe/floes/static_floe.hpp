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
                    m_thickness{1}, m_C_w{5 * 1e-3}, m_area{-1}, m_moment_cst{-1}, m_fracture{}, m_ind_fracture{} {}
                    
    StaticFloe(geometry_type new_border) : m_frame{0,0,0}, m_geometry{new_border}, m_mesh{nullptr}, m_density{917}, m_mu_static{0.7},
                    m_thickness{1}, m_C_w{5 * 1e-3}, m_area{-1}, m_moment_cst{-1}, m_fracture{}, m_ind_fracture{} {}                    
                   

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

    //! fracture accessor
    inline std::vector<point_type> get_fracture() const { return m_fracture; }
    inline void set_fracture(std::vector<point_type> frac){ m_fracture=frac; }
    inline std::vector<std::size_t> get_ind_fracture() const { return m_ind_fracture; }
    inline void set_ind_fracture(std::vector<std::size_t> ind_frac){ m_ind_fracture=ind_frac; }
    //inline bool is_fractured() const { return m_is_fractured; }
    //inline void set_is_fractured(bool fracture) { m_is_fractured=fracture; }

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
    
	// inline std::vector<StaticFloe> fracture_floe(); //{point_type mass_center {0,0};}
	inline void update_caracteristic(real_type init_density,real_type init_mu_static, real_type init_thickness, real_type init_C_w);
	// inline void generate_mesh(); 
	inline point_type get_mass_center() ;
	inline std::vector<TGeometry> fracture_floe();
    inline std::vector<TGeometry> fracture_floe_2();
    inline std::vector<TGeometry> fracture_floe_3();

    // tools for calculate energy
    inline real_type calculate_energy(real_type condition_dirichlet);
    inline real_type calculate_energy_fractured_floe(real_type condition_dirichlet,std::vector<size_t> triangle_fracture,point_type crack_start,point_type crack_stop);
    bool find_fracture(real_type condition_dirichlet,size_t edge_contact);


private:

    frame_type m_frame;         //!< Frame
    Uptr_geometry_type m_geometry;  //!< Geometry (border)
    mesh_type* m_mesh;          //!< Mesh
    density_type m_density;     //!< Density
    real_type m_mu_static;     //!< Static friction coefficient
    real_type m_thickness;     //!< Vertical thickness (constant over floe surface)
    real_type m_C_w;     //!< Oceanic skin drag coeff
    //bool m_is_fractured;     //! (if floe has been fractured)
    std::vector<point_type> m_fracture;
    std::vector<std::size_t> m_ind_fracture;

    mutable real_type m_area;  //!< Area (cached)
    // mutable density_type m_mass; //!< Mass (cached)
    mutable real_type m_moment_cst; //!< Momentum constant (cached)

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


};

template <typename T,typename TPoint,typename TGeometry,typename TMesh,typename TFrame ,typename TDensity>
std::vector<TGeometry> 
StaticFloe<T,TPoint,TGeometry,TMesh,TFrame,TDensity>::fracture_floe()
{
	// arbitrary choice, there is 24 edges so I choose to have the fracture
    // crossed on the opposite edge
    std::size_t crack_it = std::size_t(this->geometry().outer().size() / 2);
    point_type border1_begin {this->geometry().outer()[0] + (this->geometry().outer()[1] - this->geometry().outer()[0]) * 0.6};
    point_type border1_end {this->geometry().outer()[crack_it] + (this->geometry().outer()[crack_it + 1] - this->geometry().outer()[crack_it]) * 0.4};
    point_type border2_begin {this->geometry().outer()[crack_it] + (this->geometry().outer()[crack_it + 1] - this->geometry().outer()[crack_it]) * 0.6};
    point_type border2_end {this->geometry().outer()[0] + (this->geometry().outer()[1] - this->geometry().outer()[0]) * 0.4};
	// create geometry type...
	geometry_type border1;
	// rajouter des point entre les deux !
	border1.outer().push_back(border1_begin);
	for (std::size_t i = 1; i <= crack_it ; ++i){
		border1.outer().push_back(this->geometry().outer()[i]);
	}
	border1.outer().push_back(border1_end);
	
	geometry_type border2;
	border2.outer().push_back(border2_begin);
	for (std::size_t i = crack_it + 1; i < this->geometry().outer().size() ; ++i){
		border2.outer().push_back(this->geometry().outer()[i]);
	}
    border2.outer().push_back(this->geometry().outer()[0]);
	border2.outer().push_back(border2_end);
	
	// create two floe, generate by the border and caracteristic of the initial floe ?
	std::vector<TGeometry> new_border {border1,border2};
    return new_border;
}


template <typename T,typename TPoint,typename TGeometry,typename TMesh,typename TFrame ,typename TDensity>
std::vector<TGeometry> 
StaticFloe<T,TPoint,TGeometry,TMesh,TFrame,TDensity>::fracture_floe_2()
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
    real_type crack_width = 1.0; //std::sqrt(this->area()) * 0.02;
    point_type crack_ortho = direct_orthogonal(crack_start) / norm2(crack_start);
    point_type crack_delta = crack_ortho * crack_width / 2;
    crack.outer().push_back(crack_start * 2 + crack_delta);
    crack.outer().push_back(crack_start * 2 - crack_delta);
    crack.outer().push_back(- crack_start * 10 - crack_delta);
    crack.outer().push_back(- crack_start * 10 + crack_delta);
    boost::geometry::correct(crack);
    // remove crack from floe geometry
    boost::geometry::difference(this->geometry().outer(), crack, new_borders);
    return new_borders;
}


template <typename T,typename TPoint,typename TGeometry,typename TMesh,typename TFrame ,typename TDensity>
std::vector<TGeometry> 
StaticFloe<T,TPoint,TGeometry,TMesh,TFrame,TDensity>::fracture_floe_3()
{
    // create two new floe
    auto& boundary = this->geometry().outer();
    std::vector<TGeometry> new_borders;
    // crack is a long and thin rectangle containing crack_start and crack stop
    geometry_type crack;

    real_type crack_width = 0.1;//std::sqrt(this->area()) * 0.02;
    point_type crack_ortho;
    point_type crack_delta;

    // fractire [i] is on the border boundary[ind_fracture]

    for (std::size_t i = 0; i< m_ind_fracture.size() ; ++i) {
        if (m_ind_fracture[i]==boundary.size()-1){
            crack_ortho = (boundary[0]-boundary[m_ind_fracture[i]]);
            crack_ortho=crack_ortho / norm2(crack_ortho);
        }
        else {
            crack_ortho = (boundary[m_ind_fracture[i]+1]-boundary[m_ind_fracture[i]]); // / norm2(boundary[m_ind_fracture[i]+1]-boundary[m_ind_fracture[i]]);
            crack_ortho=crack_ortho / norm2(crack_ortho);
        }
        crack_delta = crack_ortho * crack_width / 2;  
        crack.outer().push_back(m_fracture[i] + crack_delta);
        crack.outer().push_back(m_fracture[i] - crack_delta);
    }
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
        
        //TPoint mass_center {0,0};
        
        return mass_center;
}


template <typename T,typename TPoint,typename TGeometry,typename TMesh,typename TFrame ,typename TDensity>
bool
StaticFloe<T,TPoint,TGeometry,TMesh,TFrame,TDensity>::find_fracture(real_type condition_dirichlet,size_t edge_contact)
{
    // first we suppose we suppose we have only 1 edge_contact

    auto& boundary = this->geometry().outer(); 
    int nb_border_edge=this->geometry().outer().size();

    // initialize all border condition as zero neumann condition
    std::vector<real_type> condition_bord(nb_border_edge,0.0); // initialize all border condition as zero neumann condition
    // get dirichlet conditions around the contact point
    condition_bord[edge_contact]=condition_dirichlet;
    
    // initial energy of the floe is energy without fracture
    real_type floe_energy {calculate_energy(condition_dirichlet)};
    real_type floe_fractured_energy {floe_energy};

    // searche all possible fracture
    point_type crack_start; // point fracture start
    point_type crack_stop; // point fracture stop
    std::size_t nb_crack_start=10; // nb of discretisation to initialise fracture
    std::size_t nb_crack_stop=10; // nb of discretisation to end fracture on each edge
    //std::vector<point_type> minimal_fract {crack_start,crack_start};

    bool is_fractured=false;

    // initialize crack starting point on the dirichlet edge
    for (std::size_t i = 1; i < nb_crack_start; ++i){  
    
        if (edge_contact<nb_border_edge-1) {
            crack_start = boundary[edge_contact]+(boundary[edge_contact+1]-boundary[edge_contact])*i/nb_crack_start;
        }
        else {// edge_contact=boundary.size()-1;
            crack_start = boundary[edge_contact]+(boundary[0]-boundary[edge_contact])*i/nb_crack_start;
        }
        // initialize crack stoping on a neumann edge
        for (std::size_t edge_stop = 0; edge_stop < nb_border_edge ; ++edge_stop){
        
            if (edge_stop != edge_contact) {
                for (std::size_t j = 1; j < nb_crack_stop; ++j){
                    if (edge_stop<nb_border_edge-1){
                        crack_stop= boundary[edge_stop]+(boundary[edge_stop+1]-boundary[edge_stop])*j/nb_crack_stop;
            // I have a fracture crack_start -> crack_stop
                    }
                    else {
                        crack_stop= boundary[edge_stop]+(boundary[0]-boundary[edge_stop])*j/nb_crack_stop;
                    }
                    /// find fractured triangle
                    std::vector<std::size_t> triangle_fract;
                    for (std::size_t k=0;k< this->m_mesh->cells().size();++k){
                        if (this->m_mesh->cells()[k].is_fractured(crack_start,crack_stop)) {
                        triangle_fract.push_back(k);
                        }
                    }
                    // ya surement moyen de faire mieux avec boost ??
                     
                    // calcul fractured floe energy
                    floe_fractured_energy=calculate_energy_fractured_floe(condition_dirichlet,triangle_fract,crack_start,crack_stop);

                        // if this energy is small, kinp in mind the fracture
                    if (floe_fractured_energy< floe_energy) {
                        floe_energy=floe_fractured_energy;
                        this->set_fracture({crack_start,crack_stop});
                        this->set_ind_fracture({edge_contact,edge_stop});
                        is_fractured=true;
                    }
                }  
            }
        }
    }
    // At this step, we check all the fracture, if there existe a fracture that minimise the total energy
    // then minimal_fracture give to {crack_start,crack_stop};
    
    return is_fractured;
}




template <typename T,typename TPoint,typename TGeometry,typename TMesh,typename TFrame ,typename TDensity>
typename StaticFloe<T,TPoint,TGeometry,TMesh,TFrame,TDensity>::real_type
StaticFloe<T,TPoint,TGeometry,TMesh,TFrame,TDensity>::calculate_energy(real_type condition_bord){ // + edge_contact

    real_type floe_energy {this->area()};

    return floe_energy;
}


template <typename T,typename TPoint,typename TGeometry,typename TMesh,typename TFrame ,typename TDensity>
typename StaticFloe<T,TPoint,TGeometry,TMesh,TFrame,TDensity>::real_type
StaticFloe<T,TPoint,TGeometry,TMesh,TFrame,TDensity>::calculate_energy_fractured_floe(real_type condition_dirichlet,std::vector<size_t> triangle_fracture,point_type crack_start,point_type crack_stop){
    real_type floe_energy;
    if (abs(condition_dirichlet)<1000.0){
        floe_energy=this->area();
    }
    else {

    floe_energy=this->area()/abs(0.001*condition_dirichlet);
    //std::cout<<" debut "<<std::endl;
    //std::cout<<this->area()<<std::endl;
    //std::cout<<condition_dirichlet<<std::endl;
    //std::cout<<norm2(crack_stop-crack_start)<<std::endl;
    real_type triangle_surface {0.0};

    for (std::size_t k=0;k< triangle_fracture.size();k++) {
       //triangle_surface = this->m_mesh->cells()[k].area();
       //triangle_surface=abs((Tri[1].x-Tri[0].x)*(Tri[2].y-Tri[0].y)-(Tri[2].x-Tri[0].x)*(Tri[1].y-Tri[0].y));
       floe_energy -= (this->m_mesh->cells()[triangle_fracture[k]].area())/abs(0.001*condition_dirichlet);
    }
    
    }

    return floe_energy+0.1*norm2(crack_stop-crack_start);
}




}} // namespace floe::floes


#endif // FLOE_FLOES_STATIC_FLOE_HPP

