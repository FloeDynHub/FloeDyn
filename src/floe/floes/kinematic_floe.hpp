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
// #include "floe/geometry/arithmetic/point_operators.hpp"

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
    using real_type = typename TStaticFloe::real_type;
    using point_type = typename TStaticFloe::point_type;
    using geometry_type = typename TStaticFloe::geometry_type;
    using mesh_type = typename TStaticFloe::mesh_type;
    using frame_type = typename TStaticFloe::frame_type;
    using static_floe_type = TStaticFloe;
    using state_type = TState;
    using floe_h_type = floe::floes::Floe_h<mesh_type>;
    using Uptr_geometry_type = std::unique_ptr<geometry_type>;
    using floe_interface_type = FloeInterface<TStaticFloe, TState>;
    using fem_problem_type = floe::fem::FemProblem<floe::floes::KinematicFloe<floe::floes::StaticFloe<real_type>>>;
    // using fem_problem_type = floe::fem::FemProblem;





    KinematicFloe() : m_geometry{nullptr}, m_floe{nullptr}, m_state{ {0,0}, 0, {0,0}, 0, {0,0}, true},
                    //   m_obstacle{false}, m_floe_h{}, m_total_impulse_received{0}, m_fem_problem{this} {}
                      m_obstacle{false}, m_floe_h{}, m_total_impulse_received{0}, m_fem_problem{this}, m_fem_problem_is_set{false} {}

    KinematicFloe(static_floe_type new_static_floe) : m_geometry{nullptr}, m_floe{new_static_floe}, m_state{ {0,0}, 0, {0,0}, 0, {0,0}, true},
                    //   m_obstacle{false}, m_floe_h{}, m_total_impulse_received{0} {} //comment je fais avec floe_h ????
                      m_obstacle{false}, m_floe_h{}, m_total_impulse_received{0}, m_fem_problem_is_set{false}
                    //   m_obstacle{false}, m_floe_h{}, m_total_impulse_received{0}, m_fem_problem{this}
    {}

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

    //! Impluse energy
    real_type impulse_energy() const;

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
    bool has_been_impacted() const {return m_total_current_impulse_received > 0;}
    // bool has_been_impacted() const {return m_total_impulse_received > 0;}

    //  std::vector<KinematicFloe<TStaticFloe,TState>> fracture_floe();
    void reset_detailed_impulse() const { m_detailed_impulse_received.clear(); }
    std::vector<geometry_type> fracture_floe();
    std::vector<geometry_type> fracture_floe_from_collisions();
    std::vector<geometry_type> fracture_floe_from_collisions_fem();

    void update_after_fracture(const state_type init_state,const bool init_obstacle_m,const real_type init_total_impulse_received, point_type mass_center_floe_init);

    bool is_active() const { return this->m_state.is_active(); }

    //! Ice speed at point p
    point_type ice_speed(point_type p) const {
        return m_state.speed + m_state.rot * fg::direct_orthogonal(p - m_state.pos);
    }

    bool set_fem_problem();
    bool prepare_elasticity();
    bool solve_elasticity();
    inline std::vector<real_type> get_fem_solution() const {return m_fem_problem.get_solution_vector();};
    // inline std::vector<real_type> get_fem_solution() {return m_fem_problem.get_solution_vector();};
    // inline std::vector<std::vector<real_type>> get_fem_stress() const {return m_fem_problem.get_stress_vector();};
    inline std::vector<real_type> get_fem_stress() const {return m_fem_problem.get_stress_vector();};
    std::size_t boundary_nb_points() const {
        return this->geometry().outer().size();
    }

    void add_contact_impulse(point_type contact_point, point_type impulse, real_type time) const;

    std::vector<point_type> get_dirichlet_condition(real_type time) const;

    real_type min_diameter() const {
        return this->static_floe().min_diameter();
    }

    bool unset_fem_problem_prepared()
    {
        m_fem_problem.unset_prepared();
        m_fem_problem_is_set = false;
        return true;
    }

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
    bool m_fem_problem_is_set;
    /*! keep track of recent collisions
     *  accumulate projected impulses on floe's boundary edges for discretized time
     *  m_recent_impulse_received will keep only recent time informations
     */
    mutable std::map<real_type, std::vector<point_type>> m_detailed_impulse_received;
    mutable std::vector<point_type> m_last_impulses;
    mutable std::vector<point_type> m_dirichlet_condition;
};

template < typename TStaticFloe, typename TState >
typename KinematicFloe<TStaticFloe,TState>::real_type
KinematicFloe<TStaticFloe,TState>::impulse_energy() const {
    auto impulsive_energy = 0;
    for (auto it = m_detailed_impulse_received.begin(); it != m_detailed_impulse_received.end(); it++) {
      auto timed_impulse = it->second;
      for (auto i = 0; i < timed_impulse.size(); i++) {
        auto impulse_norm = norm2(timed_impulse[i]);
        impulsive_energy += 0.5 * impulse_norm * impulse_norm / m_floe->mass();
      }
    }
    return impulsive_energy;
}

template < typename TStaticFloe, typename TState >
void KinematicFloe<TStaticFloe,TState>::add_contact_impulse(point_type contact_point, point_type impulse, real_type time) const {
    // round time to 1e-1
    real_type t = std::round(time * 10) / 10;
    if (m_detailed_impulse_received.find(t) == m_detailed_impulse_received.end()) {
        // not found
        std::size_t nb_points = this->boundary_nb_points();
        m_detailed_impulse_received[t] = std::vector<point_type>(nb_points, {0,0});
    }
    // find closest boundary point to contact_point
    std::size_t closest_point = 0;
    real_type min_dist = norm2(contact_point - this->geometry().outer()[0]);
    for (std::size_t i = 1; i < this->boundary_nb_points(); ++i) {
        real_type dist = norm2(contact_point - this->geometry().outer()[i]);
        if (dist < min_dist) {
            min_dist = dist;
            closest_point = i;
        }
    }
    // Add impulses at closest_point to m_detailed_impulse_received[t]
    m_detailed_impulse_received[t][closest_point] += impulse;
    add_current_impulse(norm2(impulse));

    // Remove old entries from m_detailed_impulse_received (keep only 1s)
    // (std::map is ordered ...)
    auto it = m_detailed_impulse_received.begin();
    while (it != m_detailed_impulse_received.end() && it->first < t - 1) {
        it = m_detailed_impulse_received.erase(it);
    }
    // Update dirichlet condition
    // /!\ this is not correct, because it is only update when new impulses come in, but it's not a problem for now
    // because if the condition is not enough, the floe will not fracture and the condition could only decrease
    auto dirichlet_condition = get_dirichlet_condition(time);
    // Rotate dirichlet condition for static floe
    m_dirichlet_condition.resize(dirichlet_condition.size());
    const auto trans = geometry::frame::transformer(frame_type{ {0, 0},-state().theta });
    for (auto i = 0; i < dirichlet_condition.size(); i++) {
        geometry::transform(dirichlet_condition[i], m_dirichlet_condition[i], trans);
    }
    // Update Geometry (if any)
    geometry::transform( m_floe->get_geometry(), *m_geometry, trans );


    // updating m_last_impulses, which contains only the last impulse for each mesh node.
    // m_last_impulses will be reset in solve_elasticity()
    closest_point = 0;
    auto coordinates = this->mesh().points();
    m_last_impulses.resize(coordinates.size());
    min_dist = norm2(contact_point - coordinates[0]);
    for (std::size_t i = 1; i < coordinates.size(); ++i) {
        real_type dist = norm2(contact_point - coordinates[i]);
        if (dist < min_dist) {
            min_dist = dist;
            closest_point = i;
        }
    }
    m_last_impulses[closest_point] += impulse;
}

template < typename TStaticFloe, typename TState >
std::vector<typename TStaticFloe::point_type> KinematicFloe<TStaticFloe,TState>::get_dirichlet_condition(real_type time) const {
    auto nb_points = this->boundary_nb_points();
    std::vector<point_type> resp(nb_points, {0,0});
    // iter over m_detailed_impulse_received and cumul impulses
    for (const auto t : m_detailed_impulse_received) {
        if (t.first > time - 1) {
            for (int i = 0; i < nb_points; ++i) {
                resp[i] += t.second[i];
            }
        }
    }
    return resp;
}


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
        geometry::transform( m_floe->get_mesh(), mesh(), trans ); // TODO bad_array_new_length PB HERE
    }
    // else
    // {
    //     WHEREAMI
    //     std::cout << "Warning : the updated static floe does not have a mesh. " << std::endl;
    //     std::cout << "Disabling the associated FEM Problem. " << std::endl;
    //     m_fem_problem.disable();
    // }

    // std::cout << "KinematicFloe::update() 6" << std::endl;

}

template < typename TStaticFloe, typename TState >
typename KinematicFloe<TStaticFloe,TState>::real_type
KinematicFloe<TStaticFloe,TState>::kinetic_energy() const
{
    return 0.5 * ( mass() * geometry::dot_product( m_state.speed, m_state.speed ) + moment_cst() * m_state.rot * m_state.rot );
}


template < typename TStaticFloe, typename TState >
std::vector<typename TStaticFloe::geometry_type>
KinematicFloe<TStaticFloe,TState>::fracture_floe(){
    // fracture floe (almost arbitrary fracture for now)
    return this->static_floe().fracture_floe();
}

template < typename TStaticFloe, typename TState >
std::vector<typename TStaticFloe::geometry_type>
KinematicFloe<TStaticFloe,TState>::fracture_floe_from_collisions(){
    if (this->is_obstacle()) return {}; // no fracture for obstacles
    if (this->impulse_energy() > this->static_floe().min_crack_energy()) {
        std::cout << "Floe fractured because " << this->impulse_energy() << " >= " << this->static_floe().min_crack_energy() << std::endl;
        // return this->static_floe().fracture_floe();
        return this->static_floe().fracture_floe_from_impulses(m_dirichlet_condition);
    }
    return {};
}

/**
 * @brief Looks for cracks according to fem elastic energy minimisation
 * @details if the floe is not an obstacle and did collide during the last time step,
 * the elasticity computation is
 *  * prepared if needed, and performed.
 *  * Then, lines of fracture across the floe are build by traveling along the floe contour.
 *  * If one of them corresponds to an energy-minimizing fracture, the floe is broken.
 *
 * @return list of two floes if the floe is broken, otherwise an empty list
 */
template < typename TStaticFloe, typename TState >
std::vector<typename TStaticFloe::geometry_type>
KinematicFloe<TStaticFloe,TState>::fracture_floe_from_collisions_fem(){
    // 0 - obstacles and blocks that did not collide won't break. Floes without meshes are not considered either.
    // WHEREAMI
    if (this->is_obstacle() || m_total_current_impulse_received == 0) return {};
    if (!m_floe->has_mesh()) {WHEREAMI return {};}
    if (!has_static_floe()) {WHEREAMI return {};}
    // WHEREAMI
    if (m_fem_problem.is_disabled())
    {
        std::cout << "Floe is disabled" << std::endl;
        return {};
    }
    // WHEREAMI
    // 1 - initialisation if required
    if (!prepare_elasticity())
    {
        std::cout << "FEM initialization has failed" << std::endl;
        return {};
    }
    // WHEREAMI
    // 2 - resolution
    if (!solve_elasticity())
    {
        std::cout << "FEM resolution has failed." << std::endl;
        std::cerr << "FEM resolution has failed." << std::endl;
        return {};
    }
    // WHEREAMI
    // 3 - looking for possible crack lines, and 4 - dealing with fracture
    if (m_fem_problem.get_total_elastic_energy() > 0)
    {
        real_type energy(0);
        point_type a; // crack_start
        point_type b; // crack_end
        point_type a_translated;
        point_type b_translated;
        point_type best_a;
        point_type best_b;
        point_type mass_center (get_frame().center());
        real_type theta (get_frame().theta());

        for (size_t i = 0; i < m_geometry->outer().size(); ++i)
        {
            a = m_floe->geometry().outer()[i];
            for (size_t j = i+2; j < m_geometry->outer().size()-1; ++j)
            {
                b = m_floe->geometry().outer()[j];
                if ((b != a)||(abs(int(j)-int(i)) > 3))
                {
                    a_translated = point_type (a.x*cos(theta)-a.y*sin(theta),a.x*sin(theta)+a.y*cos(theta))+mass_center;
                    b_translated = point_type (b.x*cos(theta)-b.y*sin(theta),b.x*sin(theta)+b.y*cos(theta))+mass_center;
                    real_type e = m_fem_problem.energy_release_by_breaking_along(a_translated, b_translated);
                    if (e > energy)
                    {
                        energy = e; best_a = a; best_b = b;
                    }
                }
            }
        }

        if (energy > 0)
        {
            std::cout << std::endl << "Breaking along (" << best_a.x << ";" << best_a.y << ")" << " -- (" << best_b.x << ";" << best_b.y << ")" << std::endl;
            std::cout << "Impact summary : fracturing along " << best_a << " -- " << best_b << " ; " <<  m_fem_problem.get_impact_definition() << std::endl; // this print is used to help build a database by parsing the logs 
            return this->static_floe().fracture_floe_along(best_a, best_b);
        }
        else
        {
            std::cout << "Impact summary : no fracture ; " <<  m_fem_problem.get_impact_definition() << std::endl; // idem, this print is used to help build a database by parsing the logs 
        }
    }
    return {};
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
    m_fem_problem.addFloe(this);
}


//! Update frame, geometry and mesh with respect to the current state.
template < typename TStaticFloe, typename TState >
bool
KinematicFloe<TStaticFloe,TState>::set_fem_problem()
{
    m_fem_problem.addFloe(this);
    m_fem_problem_is_set = true;
    return true;
}

template < typename TStaticFloe, typename TState >
bool
KinematicFloe<TStaticFloe,TState>::prepare_elasticity()
{
    if (!m_fem_problem_is_set)
        set_fem_problem();
        // There seems to be an issue with the floe pointer. Resetting it at each computation is necessary to prevent from segfaults. Vraisemblablement due to floe.static_floe().attach_mesh_ptr(&floe.get_floe_h().m_static_mesh); in fracture_floes().
    if (!m_fem_problem_is_set){WHEREAMI return false;}
	return m_fem_problem.prepare();
}

template < typename TStaticFloe, typename TState >
bool
KinematicFloe<TStaticFloe,TState>::solve_elasticity()
{
    if (!m_fem_problem_is_set){WHEREAMI return false;}

    bool testCase(false); // to validate code parts, using the 2D clamped-loaded beam analytical test case
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
        // ça, c'est du gros bluff en attendant une belle expression pour la CL de contact
        real_type fact_arbitraire(1E-9); // permet de conserver l'ordre de grandeur dans le cas test du bloc circulaire de r = 100m, Ecinétique avant le choc ~ Epotentielle.
        size_t nb_points = this->mesh().points().size();
        for (size_t iPoint = 0; iPoint < m_last_impulses.size(); iPoint++)
        {
            if (norm2(m_last_impulses[iPoint]) > 0)
            {
                dirichletValues.push_back(m_last_impulses[iPoint]*fact_arbitraire);
                dirichletPoints.push_back(iPoint);
            }
        }
        // resetting before next time step
        m_last_impulses.clear();
        m_last_impulses.resize(nb_points);
        if (dirichletPoints.size() > 0)
        {
            // std::cout << "Performing resolution with " << dirichletPoints.size() << " contact points" << std::endl;
            m_fem_problem.performComputation(dirichletPoints, dirichletValues);
        }
    }

    return true;
}



}} // namespace floe::floes
#endif // FLOE_FLOES_KINEMATIC_FLOE_HPP
