/*!
 * \file floe/floes/kinematic_floe.hpp
 * \brief Definition and manipulation of floe at kinematic level.
 * \author Roland Denis, Quentin Jouet
 */

#ifndef FLOE_FLOES_KINEMATIC_FLOE_HPP
#define FLOE_FLOES_KINEMATIC_FLOE_HPP

#define WHEREAMI std::cout << std::endl << "no crash until line " << __LINE__ << " in the file " __FILE__ << std::endl;



#include "floe/floes/static_floe.hpp"
#include "floe/floes/floe_exception.hpp"

#include "floe/state/space_time_state.hpp"
#include "floe/geometry/frame/frame_transformers.hpp"

#include "floe/floes/floe_h.hpp"
#include "floe/geometry/arithmetic/dot_product.hpp"
#include "floe/geometry/arithmetic/arithmetic.hpp"

#include "floe/floes/floe_interface.hpp"
#include "floe/fem/fem_problem.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>


namespace floe { namespace floes
{

/*! Kinematic of a floe
 *
 * It depends on the static of a floe.
 * Its properties (frame, geometry, mesh) are defined with respect to the absolute frame and are updated at each state change.
 *
 * Its movement is observable.
 *
 * \tparam TStaticFloe  Type of the underlying static floe.
 * \tparam TState       Type of space-time state.
 *
 */

namespace fg = floe::geometry;

template <
    typename TStaticFloe,
    typename TState = state::SpaceTimeState< typename TStaticFloe::point_type, typename TStaticFloe::real_type > 
>
class KinematicFloe : public FloeInterface<
    TStaticFloe,
    TState
>
{

public:
    
    // Type traits
    typedef typename TStaticFloe::real_type     real_type;
    typedef typename TStaticFloe::point_type    point_type;
    typedef typename TStaticFloe::geometry_type geometry_type;
    typedef typename TStaticFloe::mesh_type     mesh_type;
    typedef typename TStaticFloe::frame_type    frame_type;
    typedef TStaticFloe                         static_floe_type;
    typedef TState                              state_type;

    using floe_h_type = floe::floes::Floe_h<mesh_type>;
    using Uptr_geometry_type = std::unique_ptr<geometry_type>;
    using floe_interface_type = FloeInterface<TStaticFloe, TState>;
    using fem_problem_type = floe::fem::FemProblem<floe::floes::KinematicFloe<floe::floes::StaticFloe<real_type>>>;
    // using fem_problem_type = floe::fem::FemProblem;





    KinematicFloe() : m_geometry{nullptr}, m_floe{nullptr}, m_state{ {0,0}, 0, {0,0}, 0, {0,0}, true},
                      m_obstacle{false}, m_floe_h{}, m_total_impulse_received{0}, m_fem_problem{this} {}
                      
    KinematicFloe(static_floe_type new_static_floe) : m_geometry{nullptr}, m_floe{new_static_floe}, m_state{ {0,0}, 0, {0,0}, 0, {0,0}, true},
                    //   m_obstacle{false}, m_floe_h{}, m_total_impulse_received{0} {} //comment je fais avec floe_h ????
                      m_obstacle{false}, m_floe_h{}, m_total_impulse_received{0}, m_fem_problem{this} {} //comment je fais avec floe_h ????

    //! Deleted copy constructor
    KinematicFloe( KinematicFloe<TStaticFloe,TState> const& ) = delete;

    //! move constructor
    KinematicFloe( KinematicFloe<TStaticFloe, TState>&& ) = default;

    //! Deleted copy operator
    KinematicFloe& operator= (KinematicFloe<TStaticFloe,TState> const& ) = delete;

    //! Update geometry and mesh in respect with his current state
    void update();

    //! Static floe accessors
    inline  void                    attach_static_floe_ptr( std::unique_ptr<static_floe_type> floe )   { m_floe = std::move(floe); update(); }
    inline  static_floe_type const& static_floe()                       const   { return *m_floe; }
    inline  static_floe_type &      static_floe()                               { return *m_floe; } //!< \warning need to manually call update() after modifications.
    inline  bool                    has_static_floe()                   const   { return m_floe != nullptr; }
    inline  static_floe_type const& get_static_floe()                   const   { return *m_floe; }
    inline  void                    set_static_floe( static_floe_type const& floe )
    {
        if (! has_static_floe()) m_floe = new static_floe_type();
        *m_floe = floe;
        update();
    }

    //! Validity check (ie. has static floe)
    inline  bool    is_valid()  const { return has_static_floe(); }

    //! Frame accessors
    inline frame_type const&    get_frame()     const   { return m_floe->get_frame(); }
    inline frame_type const&    frame()         const   { return m_floe->frame(); }

    //! Geometry accessors
    inline geometry_type const&     get_geometry()      const   { return *m_geometry; }
    inline geometry_type const&     geometry()          const   { return *m_geometry; }
    inline bool                     has_geometry()      const   { return m_geometry != nullptr; }

    //! Mesh accessors
    inline mesh_type const&     get_mesh()      const   { return m_floe_h.m_kinematic_mesh; }
    inline mesh_type const&     mesh()          const   { return m_floe_h.m_kinematic_mesh; }
    inline mesh_type&     mesh()            { return m_floe_h.m_kinematic_mesh; }
    // inline bool                 has_mesh()      const   { return m_mesh != nullptr; }

    //! Area
    inline real_type area() const { return has_static_floe() ? ( m_floe->area() ) : -1; }

    //! Mass
    inline real_type mass() const { return has_static_floe() ? ( m_floe->mass() ) : -1; }

    //! State accessors
    inline state_type const&    get_state()     const   { return m_state; }
    inline void                 set_state( state_type const& state )    { 
        m_state = state;
        update();
    }
    inline state_type&    state()         const   { return m_state; }
    inline state_type &         state()                 { return m_state; } //!< \warning needs to manually call update() after modification.

    //! Obstacle
    inline bool const&  is_obstacle()   const   { return m_obstacle; }
    inline bool &       is_obstacle()           { return m_obstacle; }

    //! Momentum constant
    inline real_type moment_cst()  const { return has_static_floe() ? ( m_floe->moment_cst() ) : -1; } // Better throw an exception ...

    //! Density accessors
    inline real_type   get_density() const { return has_static_floe() ? (m_floe->get_density()) : -1; } // Throw an exception !!!

    //! Mu accessors
    inline real_type    mu_static() const   { return has_static_floe() ? m_floe->mu_static() : -1; }
    inline void set_mu_static(real_type mu_static) { if (has_static_floe()) m_floe->set_mu_static(mu_static); }

    //! Floe_h accessor
    inline floe_h_type& get_floe_h() { return m_floe_h; }

    //! Kinetic energy
    real_type kinetic_energy() const;

    //! Get total received impulse
    real_type total_received_impulse() const { return m_total_impulse_received; }
    //! Add received impulse
    void add_impulse(real_type impulse) const { m_total_impulse_received += impulse; }
    //! Reset received impulse
    void reset_impulse(real_type new_impulse = 0) const { m_total_impulse_received = new_impulse; }
    //! Idem, but recalls only the impulse received during the last impact time step 
    void add_current_impulse(real_type impulse) const { m_total_current_impulse_received += impulse; }
    //! Reset received impulse
    void reset_current_impulse(real_type new_impulse = 0) const { m_total_current_impulse_received = new_impulse; }


    //  std::vector<KinematicFloe<TStaticFloe,TState>> fracture_floe();
    std::vector<geometry_type> fracture_floe();
    
    void update_after_fracture(const state_type init_state,const bool init_obstacle_m,const real_type init_total_impulse_received, point_type mass_center_floe_init);
    
    bool is_active() const { return this->m_state.is_active(); }

    //! Ice speed at point p
    point_type ice_speed(point_type p) const {
        return m_state.speed + m_state.rot * fg::direct_orthogonal(p - m_state.pos);
    }

    bool prepare_elasticity(); 
    bool solve_elasticity(); 
    inline std::vector<real_type> get_fem_solution() const {return m_fem_problem.get_solution_vector();};
    // inline std::vector<std::vector<real_type>> get_fem_stress() const {return m_fem_problem.get_stress_vector();};
    inline std::vector<real_type> get_fem_stress() const {return m_fem_problem.get_stress_vector();};

private:

    Uptr_geometry_type m_geometry;  //!< Geometry (border)
    // mesh_type* m_mesh;          //!< Mesh
    std::unique_ptr<static_floe_type> m_floe;   //!< Static floe
    mutable state_type  m_state;        //!< Time-Space state
    bool m_obstacle;            //!< true if this floe is an obstacle

    floe_h_type m_floe_h; //!< Discretisation of the Floe
    mutable real_type m_total_impulse_received; //!< Sum all collision impulses this floe received
    mutable real_type m_total_current_impulse_received; //!< Sum all collision impulses this floe received

    fem_problem_type m_fem_problem;
};


//! Update frame, geometry and mesh with respect to the current state.
template < typename TStaticFloe, typename TState >
void
KinematicFloe<TStaticFloe,TState>::update()
{
    if (! is_valid() ) { throw( FloeException("No static floe associated.") ); }
    
    // New frame
    m_floe->set_frame( { m_state.pos, m_state.theta } );
    
    // Transformation from this new frame to absolute frame
    const auto trans = geometry::frame::transformer( m_floe->get_frame() ); 

    // Update Geometry (if any)
    if ( m_floe->has_geometry() )
    {
        if ( ! has_geometry() ) { m_geometry.reset(new geometry_type()); }
        geometry::transform( m_floe->get_geometry(), *m_geometry, trans );
    }

    // Update Mesh (if any)
    if ( m_floe->has_mesh() )
    {
        geometry::transform( m_floe->get_mesh(), mesh(), trans );
    }
}

template < typename TStaticFloe, typename TState >
typename KinematicFloe<TStaticFloe,TState>::real_type
KinematicFloe<TStaticFloe,TState>::kinetic_energy() const
{
    return 0.5 * ( mass() * geometry::dot_product( m_state.speed, m_state.speed ) + moment_cst() * m_state.rot * m_state.rot );
}


// template < typename TStaticFloe, typename TState >
// std::vector<KinematicFloe<TStaticFloe,TState>>
// KinematicFloe<TStaticFloe,TState>::fracture_floe(){
// 	// fracture floe, today arbitrary fracture
// 	std::vector<TStaticFloe> new_static_floes {this->static_floe().fracture_floe()};
// 	// create new kinematic floe from static floe
// 	std::vector<KinematicFloe<TStaticFloe,TState>>  new_floes; // {KinematicFloe<TStaticFloe,TState>(new_static_floes[0]),KinematicFloe<TStaticFloe,TState>(new_static_floes[1])};
// 	//KinematicFloe<TStaticFloe,TState> new_floes;
// 	// update 
// 	//point_type mass_center_floe_init {this->static_floe().get_mass_center()};
// 	//new_floes[0].update_after_fracture(m_state,m_obstacle,m_total_impulse_received,mass_center_floe_init);
// 	//new_floes[1].update_after_fracture(m_state,m_obstacle,m_total_impulse_received,mass_center_floe_init);
// 	//this->m_state.desactivate();
// 	return new_floes;
// }

template < typename TStaticFloe, typename TState >
std::vector<typename TStaticFloe::geometry_type>
KinematicFloe<TStaticFloe,TState>::fracture_floe(){
	// fracture floe (arbitrary fracture for now)
	return this->static_floe().fracture_floe_2();
}

//! Update frame, geometry and mesh with respect to the current state.
template < typename TStaticFloe, typename TState >
void
KinematicFloe<TStaticFloe,TState>::update_after_fracture(const state_type init_state,const bool init_obstacle,const real_type init_total_impulse_received, point_type mass_center_floe_init)
{
	TState new_state {init_state};
	new_state.pos()=init_state.pos()+this->get_mass_center-mass_center_floe_init;
	this->set_state(new_state); // the update is inside;
	this->set_obstacle(init_obstacle);
	this->set_total_impulse_received(init_total_impulse_received);
	this->update(); // update frame ect.. 
}




template < typename TStaticFloe, typename TState >
bool
KinematicFloe<TStaticFloe,TState>::prepare_elasticity()
{
	m_fem_problem.prepare();
    return true; 
}
template < typename TStaticFloe, typename TState >
bool
KinematicFloe<TStaticFloe,TState>::solve_elasticity()
{
    bool testCase(false); // to validate code parts using the 2D clamped-loaded beam test case 
    if (testCase)
    {
        // // deuxPtitsRectangles : à gauche c'est 2 3 41
        // // deuxPtitsRectangles : à droite c'est 0 1 49
        // // std::vector<size_t> dirichletPoints = {2, 3, 41, 0, 1, 49};
        // // std::vector<size_t> dirichletPoints = {2, 3, 41, 0, 1, 49, 63};
        // // std::vector<size_t> dirichletPoints = {0};
        // // std::vector<size_t> dirichletPoints = {2, 122, 41, 114, 3};
        // // std::vector<size_t> dirichletPoints = {2,3,1493,492,1774,135,1773,1777,491,367,1778,41,1779,490,1051,1295,114}; // avec nb_cells = 2000 pour faire l'encastrement à gauche de la poutre 2D 
        // // std::vector<size_t> dirichletPoints = {2,3,41,55}; // avec nb_cells = 50
        // // std::vector<size_t> dirichletPoints = {2,3,41,135,114}; // avec nb_cells = 100, 200
        // // std::vector<size_t> dirichletPoints = {2,3,492,135,491,41,490,1051,114,367,1295}; // avec nb_cells = 1000
        // // std::vector<size_t> dirichletPoints = {11,29,5,20,14}; // avec nb_cells = 1000, c'est 5 noeuds d'affilée au millieu 
        // std::vector<size_t> dirichletPoints = {29,5,20}; // avec nb_cells = 1000, c'est 3 noeuds d'affilée 
        // point_type a(0,0);
        // point_type b(1,1);
        // point_type c(0,1);
        // point_type d(0,2);
        // point_type e(0,3);
        // point_type f(1,2);
        // // std::vector<point_type> dirichletValues = {b, a};
        // // std::vector<point_type> dirichletValues = {a, a, a, b, b, b};
        // // std::vector<point_type> dirichletValues = {a, a, a, a, a, a, b};
        // // std::vector<point_type> dirichletValues = {a, a, a, a, a};
        // // std::vector<point_type> dirichletValues = {a, a, a, a, a,a,a,a,a,a,a,a,a,a,a,a,a};
        // // std::vector<point_type> dirichletValues = {a, a, a, a, a,a,a,a,a,a,a};
        // // std::vector<point_type> dirichletValues = {b, b, b, b, b,b,b,b,b,b,b};
        // // std::vector<point_type> dirichletValues = {c,d,e,d,c};
        // std::vector<point_type> dirichletValues = {c,d,c};
        point_type a(0,0); 
        std::vector<size_t> dirichletPoints = {2, 3, 41}; 
        std::vector<point_type> dirichletValues = {a, a, a};
        m_fem_problem.performComputation(dirichletPoints, dirichletValues);
    }
    else 
    {
        std::vector<size_t> dirichletPoints = {}; 
        std::vector<point_type> dirichletValues = {};

        real_type fact_arbitraire(1E-6); // permet de conserver l'ordre de grandeur dans le cas test de la poutre 2D, Ecinétique dissipée en Epotentielle 
        point_type direction(-1,0);
        std::vector<size_t> contact_points = {0}; // 46 // to do : liste des noeuds en contact, à chercher dans le problem.m_priximity_detector.contact_graph ? 
        std::vector<real_type> amplitudes = {m_total_current_impulse_received*fact_arbitraire}; // to do : quelles variables ? velocity*mass a priori ? à voir avec Toai. Contact_point.relative_speed() ? Un coeff arbitraire pour conserver l'énergie cinétique du cas test ? 
        std::vector<point_type> directions = {direction}; // to do : vecteur normal au choc, vraisemblamenet un ContactPoint.frame.v()  
        size_t n_contact_points = contact_points.size();
        for (size_t iPoint = 0 ; iPoint < n_contact_points ; ++iPoint)
        {
            dirichletValues.push_back(amplitudes[iPoint]*directions[iPoint]);
            dirichletPoints.push_back(contact_points[iPoint]);
        }

        m_fem_problem.performComputation(dirichletPoints, dirichletValues);
    }
    
    return true; 
}



}} // namespace floe::floes
#endif // FLOE_FLOES_KINEMATIC_FLOE_HPP

