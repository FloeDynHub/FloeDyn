/*!
 * \file ope/dynamics_manager.hpp
 * \brief dynamics manager
 * \author Quentin Jouet
 */

#ifndef OPE_DYNAMICS_MANAGER_HPP
#define OPE_DYNAMICS_MANAGER_HPP

#include "floe/integration/gauss_legendre.hpp"
#include "floe/integration/integrate.hpp"
#include "floe/ope/external_forces.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

 #include <iostream> // DEBUG


namespace floe { namespace ope
{

/*! DynamicsManager
 *
 * Operator for dynamics processing
 *
 */


template <typename TFloeGroup>
class DynamicsManager
{

public:

    using floe_group_type = TFloeGroup;
    using floe_type = typename floe_group_type::floe_type;
    using point_type = typename floe_type::point_type;
    using value_type = typename floe_type::value_type;
    using integration_strategy = floe::integration::RefGaussLegendre<value_type,2,2>;
    using external_forces_type = ExternalForces<TFloeGroup>;

    DynamicsManager() = default;

    void move_floes(floe_group_type& floe_group, value_type delta_t);

protected:

    external_forces_type m_external_forces;

    virtual void move_floe(floe_type& floe, value_type delta_t);
    void translate_floe(floe_type& floe, value_type delta_t);
    void rotate_floe(floe_type& floe, value_type delta_t);
};


template <typename TFloeGroup>
void
DynamicsManager<TFloeGroup>::move_floes(floe_group_type& floe_group, value_type delta_t)
{   
    // OpenMP doesn't like this syntax
    // for (auto& floe : floe_group.get_floes())
    //     move_floe(floe, delta_t);
    // omp_set_num_threads(2);
    #pragma omp parallel for
    for (std::size_t i=0; i < floe_group.get_floes().size(); ++i)
        move_floe(floe_group.get_floes()[i], delta_t);
}


template <typename TFloeGroup>
void
DynamicsManager<TFloeGroup>::move_floe(floe_type& floe, value_type delta_t)
{
    translate_floe(floe, delta_t);
    rotate_floe(floe, delta_t);
    floe.update(); // alternative : pour ne pas avoir à update() : utiliser set_state()
}


template <typename TFloeGroup>
void
DynamicsManager<TFloeGroup>::translate_floe(floe_type& floe, value_type delta_t)
{
    auto drag_force = floe::integration::integrate(
        m_external_forces.total_drag(floe),
        floe.mesh(),
        integration_strategy()
    );
    floe.state().pos += delta_t * floe.state().speed;
    floe.state().speed += ( delta_t / floe.mass() ) * drag_force
                         + delta_t * m_external_forces.coriolis_effect(floe);
}


template <typename TFloeGroup>
void
DynamicsManager<TFloeGroup>::rotate_floe(floe_type& floe, value_type delta_t)
{
    auto rot_drag_force = floe::integration::integrate(
        m_external_forces.total_rot_drag(floe),
        floe.mesh(),
        integration_strategy()
    );
    floe.state().theta += delta_t * floe.state().rot;
    floe.state().rot += ( delta_t / floe.moment_cst() ) * rot_drag_force;
}


}} // namespace floe::ope


#endif // OPE_DYNAMICS_MANAGER_HPP