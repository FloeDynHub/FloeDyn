/*!
 * \file dynamics/periodic_dynamics_manager.h
 * \brief Extended dynamics manager handling translations for periodic space
 * \author Quentin Jouet
 */

#ifndef OPE_PERIODIC_DYNAMICS_MANAGER_H
#define OPE_PERIODIC_DYNAMICS_MANAGER_H

#include "floe/dynamics/dynamics_manager.h"


namespace floe { namespace dynamics
{

/*! DynamicsManager
 *
 * Operator for dynamics processing
 *
 */


template <typename TExternalForces, typename TFloeGroup, typename TSpaceTopology>
class PeriodicDynamicsManager : public DynamicsManager<TExternalForces, TFloeGroup>
{

public:
    using base_class = DynamicsManager<TExternalForces, TFloeGroup>;
    using floe_type = typename base_class::floe_type;
    using point_type = typename base_class::point_type;
    using value_type = typename base_class::value_type;
    using topology_type = TSpaceTopology;
    // using integration_strategy = typename base_class::integration_strategy;

    PeriodicDynamicsManager(value_type const& time_ref) :
        base_class(time_ref), m_topology{nullptr} {}

    //! Set topology
    inline void set_topology(topology_type const& t) { m_topology = &t; }

private:
    topology_type const* m_topology; //!< Space topology
    virtual void move_floe(floe_type& floe, value_type delta_t) override;
    //! Translate floe if needed according to periodic boundary conditions (topology)
    bool replace_floe(floe_type& floe);
    virtual value_type ocean_window_area() override { return m_topology->area(); }
};


}} // namespace floe::dynamics


#endif // OPE_PERIODIC_DYNAMICS_MANAGER_H