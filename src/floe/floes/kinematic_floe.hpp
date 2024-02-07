/*!
 * \file floe/floes/kinematic_floe.hpp
 * \brief Definition and manipulation of floe at kinematic level.
 * \author Roland Denis, Quentin Jouet
 */

#ifndef FLOE_FLOES_KINEMATIC_FLOE_HPP
#define FLOE_FLOES_KINEMATIC_FLOE_HPP


#include "floe/floes/static_floe.hpp"
#include "floe/floes/floe_exception.hpp"

#include "floe/state/space_time_state.hpp"
#include "floe/geometry/frame/frame_transformers.hpp"

#include "floe/floes/floe_h.hpp"
#include "floe/geometry/arithmetic/dot_product.hpp"
#include "floe/geometry/arithmetic/arithmetic.hpp"
// #include "floe/geometry/arithmetic/point_operators.hpp"

#include "floe/floes/floe_interface.hpp"

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

    KinematicFloe() : m_geometry{nullptr}, m_floe{nullptr}, m_state{ {0,0}, 0, {0,0}, 0, {0,0}, true},
                      m_obstacle{false}, m_floe_h{}, m_total_impulse_received{0} {}

    KinematicFloe(static_floe_type new_static_floe) : m_geometry{nullptr}, m_floe{new_static_floe}, m_state{ {0,0}, 0, {0,0}, 0, {0,0}, true},
                      m_obstacle{false}, m_floe_h{}, m_total_impulse_received{0} {} //comment je fais avec floe_h ????

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

    //  std::vector<KinematicFloe<TStaticFloe,TState>> fracture_floe();
    std::vector<geometry_type> fracture_floe();

    void update_after_fracture(const state_type init_state,const bool init_obstacle_m,const real_type init_total_impulse_received, point_type mass_center_floe_init);

    bool is_active() const { return this->m_state.is_active(); }

    //! Ice speed at point p
    point_type ice_speed(point_type p) const {
        return m_state.speed + m_state.rot * fg::direct_orthogonal(p - m_state.pos);
    }

    std::size_t boundary_nb_points() const {
        return this->geometry().outer().size();
    }

    void add_contact_impulse(point_type contact_point, point_type impulse, real_type time) const;

    std::vector<point_type> get_dirichlet_condition(real_type time) const;

    real_type min_radius() const {
        return this->static_floe().min_radius();
    }

private:

    Uptr_geometry_type m_geometry;  //!< Geometry (border)
    // mesh_type* m_mesh;          //!< Mesh
    std::unique_ptr<static_floe_type> m_floe;   //!< Static floe
    mutable state_type  m_state;        //!< Time-Space state
    bool m_obstacle;            //!< true if this floe is an obstacle

    floe_h_type m_floe_h; //!< Discretisation of the Floe
    mutable real_type m_total_impulse_received; //!< Sum all collision impulses this floe received

    /*! keep track of recent collisions
     *  accumulate projected impulses on floe's boundary edges for discretized time
     *  m_recent_impulse_received will keep only recent time informations
     */
    mutable std::map<real_type, std::vector<point_type>> m_detailed_impulse_received;
};

template < typename TStaticFloe, typename TState >
typename KinematicFloe<TStaticFloe,TState>::real_type
KinematicFloe<TStaticFloe,TState>::impulse_energy() const {
    auto impulsive_energy = 0;
    for (auto it = m_detailed_impulse_received.begin(); it != m_detailed_impulse_received.end(); it++) {
      auto timed_impulse = it->second;
      for (auto i = 0; i < timed_impulse.size(); i++) {
        impulsive_energy += 0.5*norm2(timed_impulse[i])^2 / m_floe->mass();
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
    // Remove old entries from m_detailed_impulse_received (keep only 1s)
    for (auto it = m_detailed_impulse_received.begin(); it != m_detailed_impulse_received.end(); ) {
        if (it->first < t - 1) {
            it = m_detailed_impulse_received.erase(it);
        } else {
            ++it;
        }
    }
    // TODO smart method for filtering keys ?
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
        geometry::transform( m_floe->get_mesh(), mesh(), trans );
    }
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



}} // namespace floe::floes
#endif // FLOE_FLOES_KINEMATIC_FLOE_HPP

