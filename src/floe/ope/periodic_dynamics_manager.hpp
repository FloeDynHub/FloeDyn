/*!
 * \file ope/periodic_dynamics_manager.hpp
 * \brief Extended dynamics manager handling translations for periodic space
 * \author Quentin Jouet
 */

#ifndef OPE_PERIODIC_DYNAMICS_MANAGER_HPP
#define OPE_PERIODIC_DYNAMICS_MANAGER_HPP

#include "floe/ope/periodic_dynamics_manager.h"


namespace floe { namespace ope
{


template <typename TExternalForces, typename TFloeGroup, typename TSpaceTopology>
void
PeriodicDynamicsManager<TExternalForces, TFloeGroup, TSpaceTopology>::move_floe(floe_type& floe, value_type delta_t)
{
    base_class::move_floe(floe, delta_t);
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


}} // namespace floe::ope


#endif // OPE_PERIODIC_DYNAMICS_MANAGER_HPP