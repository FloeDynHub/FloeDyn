/*!
 * \file floe/floes/kinematic_floe.hpp
 * \brief Definition and manipulation of floe at kinematic level.
 * \author Roland Denis
 */

#ifndef FLOE_FLOES_KINEMATIC_FLOE_HPP
#define FLOE_FLOES_KINEMATIC_FLOE_HPP


#include "floe/floes/static_floe.hpp"
#include "floe/floes/floe_exception.hpp"

#include "floe/state/space_time_state.hpp"
#include "floe/geometry/frame/frame_transformers.hpp"

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
 * \todo use smart pointer.
 */
template <
    typename TStaticFloe,
    typename TState = state::SpaceTimeState< typename TStaticFloe::point_type, typename TStaticFloe::value_type > 
>
class KinematicFloe
{

public:
    
    // Type traits
    typedef typename TStaticFloe::value_type    value_type;
    typedef typename TStaticFloe::point_type    point_type;
    typedef typename TStaticFloe::geometry_type geometry_type;
    typedef typename TStaticFloe::mesh_type     mesh_type;
    typedef typename TStaticFloe::frame_type    frame_type;
    typedef TStaticFloe                         static_floe_type;
    typedef TState                              state_type;

    //! Default constructor
    KinematicFloe() : m_geometry{nullptr}, m_mesh{nullptr}, m_floe{nullptr}, m_state{ {0,0}, 0, {0,0}, 0 }, m_obstacle{false} {}

    //! Destructor
    ~KinematicFloe() { delete m_geometry; delete m_mesh; delete m_floe; }

    //! Deleted copy constructor
    KinematicFloe( KinematicFloe<TStaticFloe,TState> const& ) = delete;

    //! rvalue reference constructor (to allow move semantics)
    KinematicFloe( KinematicFloe<TStaticFloe, TState>&& ){};

    //! Deleted copy operator
    KinematicFloe& operator= (KinematicFloe<TStaticFloe,TState> const& ) = delete;

    //! Update geometry and mesh in respect with his current state
    void update();

    //! Static floe accessors
    inline  void                    attach_static_floe_ptr( static_floe_type* floe )   { delete m_floe; m_floe = floe; update(); }
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
    inline mesh_type const&     get_mesh()      const   { return *m_mesh; }
    inline mesh_type const&     mesh()          const   { return *m_mesh; }
    inline bool                 has_mesh()      const   { return m_mesh != nullptr; }

    //! Area
    inline value_type area() const { return has_static_floe() ? ( m_floe->area() ) : -1; }

    //! Mass
    inline value_type mass() const { return has_static_floe() ? ( m_floe->mass() ) : -1; }

    //! State accessors
    inline state_type const&    get_state()     const   { return m_state; }
    inline void                 set_state( state_type const& state )    { m_state = state; update(); }
    inline state_type const&    state()         const   { return m_state; }
    inline state_type &         state()                 { return m_state; } //!< \warning needs to manually call update() after modification.

    //! Obstacle
    inline bool const&  is_obstacle()   const   { return m_obstacle; }
    inline bool &       is_obstacle()           { return m_obstacle; }

    //! Momentum constant
    inline value_type moment_cst()  const { return has_static_floe() ? ( m_floe->moment_cst() ) : -1; } // Better throw an exception ...

    //! Density accessors
    inline value_type   get_density() const { return has_static_floe() ? (m_floe->get_density()) : -1; } // Throw an exception !!!

    //! Mu accessors
    inline value_type    mu_static() const   { return has_static_floe() ? m_floe->mu_static() : -1; }

private:

    geometry_type* m_geometry;  //!< Geometry (border)
    mesh_type* m_mesh;          //!< Mesh
    static_floe_type* m_floe;   //!< Static floe
    state_type  m_state;        //!< Time-Space state
    bool m_obstacle;            //!< true if this floe is an obstacle

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
        if ( ! has_geometry() ) { m_geometry = new geometry_type(); }
        geometry::transform( m_floe->get_geometry(), *m_geometry, trans );
    }

    // Update Mesh (if any)
    if ( m_floe->has_mesh() )
    {
        if ( ! has_mesh() ) { m_mesh = new mesh_type(); }
        geometry::transform( m_floe->get_mesh(), *m_mesh, trans );
    }
}

}} // namespace floe::floes
#endif // FLOE_FLOES_KINEMATIC_FLOE_HPP

