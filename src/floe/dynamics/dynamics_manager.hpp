/*!
 * \file dynamics/dynamics_manager.hpp
 * \brief Dynamics manager
 * \author Quentin Jouet
 */

#ifndef OPE_DYNAMICS_MANAGER_HPP
#define OPE_DYNAMICS_MANAGER_HPP

#include "floe/dynamics/dynamics_manager.h"

#include "floe/integration/gauss_legendre.hpp"
#include "floe/integration/integrate.hpp"
#include "floe/io/matlab/pze_import.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

 #include <iostream> // DEBUG


namespace floe { namespace dynamics
{

template<typename T>
using integration_strategy = floe::integration::RefGaussLegendre<T,2,2>;

template <typename TExternalForces, typename TFloeGroup>
void
DynamicsManager<TExternalForces, TFloeGroup>::move_floes(floe_group_type& floe_group, value_type delta_t)
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


template <typename TExternalForces, typename TFloeGroup>
void
DynamicsManager<TExternalForces, TFloeGroup>::move_floe(floe_type& floe, value_type delta_t)
{
    if (floe.is_obstacle()) return; // Obstacles don't move

    state_type new_state = floe.state();

    // Translation part
    auto drag_force = floe::integration::integrate(
        m_external_forces.total_drag(floe),
        floe.mesh(),
        integration_strategy<value_type>()
    );
    new_state.pos += delta_t * floe.state().speed;
    new_state.speed += ( delta_t / floe.mass() ) * drag_force
                         + delta_t * m_external_forces.coriolis_effect(floe);

    // Rotation part
    auto rot_drag_force = floe::integration::integrate(
        m_external_forces.total_rot_drag(floe),
        floe.mesh(),
        integration_strategy<value_type>()
    );
    new_state.theta += delta_t * floe.state().rot;
    new_state.rot += ( delta_t / floe.moment_cst() ) * rot_drag_force;

    /* Adding random perturbation to speed and rot
       (improve collision computing, physically justifiable) */
    value_type rand_norm = 1e-6 * delta_t;
    auto dist_rot = std::uniform_real_distribution<value_type>{-rand_norm, rand_norm};
    auto dist_angle = std::uniform_real_distribution<value_type>{0, 2 * M_PI};
    auto rand_theta = dist_angle(this->m_random_generator);
    auto rand_rot = dist_rot(this->m_random_generator);
    auto rand_speed = rand_norm * point_type{cos(rand_theta), sin(rand_theta)};
    new_state.speed += rand_speed;
    new_state.rot += rand_rot;

    // Floe update
    floe.set_state(new_state);
}


template <typename TExternalForces, typename TFloeGroup>
void
DynamicsManager<TExternalForces, TFloeGroup>::update_ocean(floe_group_type& floe_group, value_type delta_t)
{   
    point_type diff_speed{0,0};
    if (m_OBL_status)
    {
        value_type floes_area = floe_group.total_area();
        value_type win_area = ocean_window_area();
        value_type water_area = win_area - floes_area;
        point_type floe_group_mass_center = floe_group.mass_center();
        value_type OBL_mass = win_area * m_external_forces.OBL_surface_mass();
        auto strategy = integration_strategy<value_type>();
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

template <typename TExternalForces, typename TFloeGroup>
void
DynamicsManager<TExternalForces, TFloeGroup>::load_matlab_ocean_window_data(std::string const& filename, floe_group_type const& floe_group)
{
    m_ocean_window_area = floe::io::matlab::ocean_window_area_from_file(filename);
    if (m_ocean_window_area == 0)
    {
        m_ocean_window_area = floe_group.bounding_window_area();
    }
}


}} // namespace floe::dynamics


#endif // OPE_DYNAMICS_MANAGER_HPP