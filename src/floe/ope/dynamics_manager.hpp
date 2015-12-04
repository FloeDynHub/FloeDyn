/*!
 * \file ope/dynamics_manager.hpp
 * \brief Dynamics manager
 * \author Quentin Jouet
 */

#ifndef OPE_DYNAMICS_MANAGER_HPP
#define OPE_DYNAMICS_MANAGER_HPP

#include "floe/integration/gauss_legendre.hpp"
#include "floe/integration/integrate.hpp"
#include "floe/io/matlab/pze_import.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

 #include <iostream> // DEBUG


namespace floe { namespace ope
{

/*! DynamicsManager
 *
 * Operator for dynamics processing (Floes and ocean kinematics)
 *
 */


template <typename TExternalForces>
class DynamicsManager
{

public:
    using external_forces_type = TExternalForces;
    using floe_group_type = typename external_forces_type::floe_group_type;
    using floe_type = typename floe_group_type::floe_type;
    using point_type = typename floe_type::point_type;
    using value_type = typename floe_type::value_type;
    using state_type = typename floe_type::state_type;
    using integration_strategy = floe::integration::RefGaussLegendre<value_type,2,2>;

    //! Constructor
    DynamicsManager(value_type const& time_ref) : m_external_forces{time_ref}, m_ocean_window_area{0}, m_OBL_status{0} {}

    //! Floes state update
    void move_floes(floe_group_type& floe_group, value_type delta_t);
    //! Ocean state update
    void update_ocean(floe_group_type& floe_group, value_type delta_t);

    //! Load ocean and wind data from a topaz file
    inline void load_matlab_topaz_data(std::string const& filename) {
        m_external_forces.load_matlab_topaz_data(filename);
    }

    //! Load ocean window area (box surrounding floes)
    void load_matlab_ocean_window_data(std::string const& filename, floe_group_type const& floe_group);

    //! for output
    inline point_type OBL_speed() const { return m_external_forces.OBL_speed(); }
    inline void set_OBL_speed(point_type OBL_speed) { return m_external_forces.update_water_speed(OBL_speed); }
    inline void set_OBL_status(int status) { m_OBL_status = status; }

    //! Accessor for specific use
    external_forces_type& get_external_forces() { return m_external_forces; }

protected:

    external_forces_type m_external_forces; //! External forces manager
    value_type m_ocean_window_area; //! Ocean window area (for OBL computing)
    int m_OBL_status; //! OBL (Oceanic Boundary Layer) status : 0 = no coupling, 1 = coupling

    //! Move one floe
    virtual void move_floe(floe_type& floe, value_type delta_t);
    //! Ocean window area accessor
    virtual value_type ocean_window_area() { return m_ocean_window_area; }
};


template <typename TExternalForces>
void
DynamicsManager<TExternalForces>::move_floes(floe_group_type& floe_group, value_type delta_t)
{   
    // OpenMP doesn't like this syntax
    // for (auto& floe : floe_group.get_floes())
    //     move_floe(floe, delta_t);
    #pragma omp parallel for
    for (std::size_t i=0; i < floe_group.get_floes().size(); ++i){
        move_floe(floe_group.get_floes()[i], delta_t);
    }
    update_ocean(floe_group, delta_t);
}


template <typename TExternalForces>
void
DynamicsManager<TExternalForces>::move_floe(floe_type& floe, value_type delta_t)
{
    if (floe.is_obstacle()) return; // Obstacles don't move

    state_type new_state = floe.state();

    // Translation part
    auto drag_force = floe::integration::integrate(
        m_external_forces.total_drag(floe),
        floe.mesh(),
        integration_strategy()
    );
    new_state.pos += delta_t * floe.state().speed;
    new_state.speed += ( delta_t / floe.mass() ) * drag_force
                         + delta_t * m_external_forces.coriolis_effect(floe);

    // Rotation part
    auto rot_drag_force = floe::integration::integrate(
        m_external_forces.total_rot_drag(floe),
        floe.mesh(),
        integration_strategy()
    );
    new_state.theta += delta_t * floe.state().rot;
    new_state.rot += ( delta_t / floe.moment_cst() ) * rot_drag_force;    

    // Floe update
    floe.set_state(new_state);
}


template <typename TExternalForces>
void
DynamicsManager<TExternalForces>::update_ocean(floe_group_type& floe_group, value_type delta_t)
{   
    point_type diff_speed{0,0};
    if (m_OBL_status)
    {
        value_type floes_area = floe_group.total_area();
        value_type win_area = ocean_window_area();
        value_type water_area = win_area - floes_area;
        point_type floe_group_mass_center = floe_group.mass_center();
        value_type OBL_mass = win_area * m_external_forces.OBL_surface_mass();
        auto strategy = integration_strategy();
        // calculate floes action on ocean
        point_type floes_force = {0, 0};
        for (auto& floe : floe_group.get_floes())
            floes_force += floe::integration::integrate(m_external_forces.ocean_drag_2(floe), floe.mesh(), strategy);
        // calculate water speed delta
        diff_speed = delta_t * ( 
            ( 1 / OBL_mass ) * ( - floes_force + water_area * m_external_forces.air_drag_ocean() )
            + m_external_forces.ocean_coriolis(floe_group_mass_center)
            + m_external_forces.deep_ocean_friction()
        );
    }
    // update water speed
    m_external_forces.update_water_speed( diff_speed );
}

template <typename TExternalForces>
void
DynamicsManager<TExternalForces>::load_matlab_ocean_window_data(std::string const& filename, floe_group_type const& floe_group)
{
    m_ocean_window_area = floe::io::matlab::ocean_window_area_from_file(filename);
    if (m_ocean_window_area == 0)
    {
        auto a = floe_group.bounding_window();
        m_ocean_window_area = (a[1] - a[0]) * (a[3] - a[2]);
    }
}


}} // namespace floe::ope


#endif // OPE_DYNAMICS_MANAGER_HPP