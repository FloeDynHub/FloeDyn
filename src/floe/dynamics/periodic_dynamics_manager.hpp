/*!
 * \file dynamics/periodic_dynamics_manager.hpp
 * \brief Extended dynamics manager handling translations for periodic space
 * \author Quentin Jouet
 */

#ifndef OPE_PERIODIC_DYNAMICS_MANAGER_HPP
#define OPE_PERIODIC_DYNAMICS_MANAGER_HPP

#include "floe/dynamics/periodic_dynamics_manager.h"


namespace floe { namespace dynamics
{


template <typename TExternalForces, typename TFloeGroup, typename TSpaceTopology>
void
PeriodicDynamicsManager<TExternalForces, TFloeGroup, TSpaceTopology>::move_floe(floe_type& floe, real_type delta_t, real_type t, bool mode_eight, bool is_first_floe)
{
    base_class::move_floe(floe, delta_t, t, mode_eight, is_first_floe);
    auto replaced = replace_floe(floe);
    if (replaced) floe.update();
}

template <typename TExternalForces, typename TFloeGroup, typename TSpaceTopology>
bool
PeriodicDynamicsManager<TExternalForces, TFloeGroup, TSpaceTopology>::replace_floe(floe_type& floe)
{
    auto trans = m_topology->replace(floe.state().pos);
    floe.state().trans += trans;
    return (trans != point_type{0,0});
}


}} // namespace floe::dynamics


#endif // OPE_PERIODIC_DYNAMICS_MANAGER_HPP