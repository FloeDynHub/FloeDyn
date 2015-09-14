/*!
 * \file ope/periodic_dynamics_manager.hpp
 * \brief dynamics manager managing translations for periodic space
 * \author Quentin Jouet
 */

#ifndef OPE_PERIODIC_DYNAMICS_MANAGER_HPP
#define OPE_PERIODIC_DYNAMICS_MANAGER_HPP

#include "floe/ope/dynamics_manager.hpp"
// #include "floe/integration/integrate.hpp"

namespace floe { namespace ope
{

/*! DynamicsManager
 *
 * Operator for dynamics processing
 *
 */


template <typename TFloeGroup, typename TSpaceTopology>
class PeriodicDynamicsManager : public DynamicsManager<TFloeGroup>
{

public:
    using base_class = DynamicsManager<TFloeGroup>;
    using floe_type = typename TFloeGroup::floe_type;
    using point_type = typename floe_type::point_type;
    using value_type = typename floe_type::value_type;
    using topology_type = TSpaceTopology;
    using integration_strategy = typename base_class::integration_strategy;

    PeriodicDynamicsManager(value_type const& time_ref) :
        base_class(time_ref), m_topology{nullptr} {}

    inline void set_topology(topology_type const& t) { m_topology = &t; }

    virtual void update_ocean(TFloeGroup& floe_group, value_type delta_t);

private:
    topology_type const* m_topology; //!< Space topology
    virtual void move_floe(floe_type& floe, value_type delta_t);
    void replace_floe(floe_type& floe);
};


template <typename TFloeGroup, typename TSpaceTopology>
void
PeriodicDynamicsManager<TFloeGroup, TSpaceTopology>::move_floe(floe_type& floe, value_type delta_t)
{
    // base_class::move_floe(floe, delta_t);
    base_class::translate_floe(floe, delta_t);
    base_class::rotate_floe(floe, delta_t);
    replace_floe(floe);
    floe.update(); // alternative : pour ne pas avoir à update() : utiliser set_state()
}

template <typename TFloeGroup, typename TSpaceTopology>
void
PeriodicDynamicsManager<TFloeGroup, TSpaceTopology>::replace_floe(floe_type& floe)
{
    m_topology->replace(floe.state().pos, floe.state().trans);
}

template <typename TFloeGroup, typename TSpaceTopology>
void
PeriodicDynamicsManager<TFloeGroup, TSpaceTopology>::update_ocean(TFloeGroup& floe_group, value_type delta_t)
{   
    value_type floes_area = floe_group.total_area();
    value_type window_area = m_topology->area();
    value_type water_area = window_area - floes_area;
    point_type window_center = m_topology->center();
    value_type OBL_mass = window_area * base_class::m_external_forces.OBL_surface_mass();
    auto strategy = integration_strategy();
    // calculate floes action on ocean
    point_type floes_force = {0, 0};
    for (auto& floe : floe_group.get_floes())
        floes_force += floe::integration::integrate(base_class::m_external_forces.ocean_drag_2(floe), floe.mesh(), strategy);
    // calculate water speed delta
    point_type diff_speed = delta_t * ( 
        ( 1 / OBL_mass ) * ( floes_force + water_area * base_class::m_external_forces.air_drag_ocean() )
        + base_class::m_external_forces.ocean_coriolis(window_center)
        + base_class::m_external_forces.deep_ocean_friction()
    );
    // update water speed
    base_class::m_external_forces.update_water_speed( diff_speed );
}



}} // namespace floe::ope


#endif // OPE_PERIODIC_DYNAMICS_MANAGER_HPP