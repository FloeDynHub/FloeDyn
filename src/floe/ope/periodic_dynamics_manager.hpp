/*!
 * \file ope/periodic_dynamics_manager.hpp
 * \brief dynamics manager managing translations for periodic space
 * \author Quentin Jouet
 */

#ifndef OPE_PERIODIC_DYNAMICS_MANAGER_HPP
#define OPE_PERIODIC_DYNAMICS_MANAGER_HPP

#include "floe/ope/dynamics_manager.hpp"

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

    PeriodicDynamicsManager(value_type const& time_ref) :
        base_class(time_ref), m_topology{nullptr} {}

    inline void set_topology(topology_type const& t) { m_topology = &t; }

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
    m_topology->replace(floe.state().pos);
}



}} // namespace floe::ope


#endif // OPE_PERIODIC_DYNAMICS_MANAGER_HPP