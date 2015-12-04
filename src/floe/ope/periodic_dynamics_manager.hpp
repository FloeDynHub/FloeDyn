/*!
 * \file ope/periodic_dynamics_manager.hpp
 * \brief Extended dynamics manager handling translations for periodic space
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


template <typename TExternalForces, typename TSpaceTopology>
class PeriodicDynamicsManager : public DynamicsManager<TExternalForces>
{

public:
    using base_class = DynamicsManager<TExternalForces>;
    using floe_type = typename base_class::floe_type;
    using point_type = typename base_class::point_type;
    using value_type = typename base_class::value_type;
    using topology_type = TSpaceTopology;
    using integration_strategy = typename base_class::integration_strategy;

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


template <typename TExternalForces, typename TSpaceTopology>
void
PeriodicDynamicsManager<TExternalForces, TSpaceTopology>::move_floe(floe_type& floe, value_type delta_t)
{
    base_class::move_floe(floe, delta_t);
    auto replaced = replace_floe(floe);
    if (replaced) floe.update();
}

template <typename TExternalForces, typename TSpaceTopology>
bool
PeriodicDynamicsManager<TExternalForces, TSpaceTopology>::replace_floe(floe_type& floe)
{
    auto trans = m_topology->replace(floe.state().pos);
    floe.state().trans += trans;
    return (trans != point_type{0,0});
}



}} // namespace floe::ope


#endif // OPE_PERIODIC_DYNAMICS_MANAGER_HPP