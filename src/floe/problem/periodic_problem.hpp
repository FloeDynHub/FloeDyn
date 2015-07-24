/*!
 * \file problem/periodic_problem.hpp
 * \brief Smooth periodic problem
 * \author Quentin Jouet
 */

#ifndef PROBLEM_PERIODIC_PROBLEM_HPP
#define PROBLEM_PERIODIC_PROBLEM_HPP

#include "floe/problem/problem.hpp"

#include <iostream> // debug

namespace floe { namespace problem
{

/*! Problem
 *
 * It represents the whole problem of moving N floes in interval time [0, T], in a space with periodic border conditions.
 *
 * \tparam TFloeGroup  
 * \tparam TProxymityDetector 
 * \tparam TCollisionManager   
 * \tparam TDynamicsManager 
 * \tparam TDomain
 * /tparam TSpaceTopology
 *
 */

template <
    typename TFloeGroup,
    typename TProxymityDetector,
    typename TCollisionManager,
    typename TDynamicsManager,
    typename TDomain,
    typename TSpaceTopology
>
class PeriodicProblem : public Problem<TFloeGroup, TProxymityDetector, TCollisionManager, TDynamicsManager, TDomain>
{
public:
    using base_class = Problem<TFloeGroup, TProxymityDetector, TCollisionManager, TDynamicsManager, TDomain>;
    using value_type = typename TFloeGroup::floe_type::value_type;
    using point_type = typename TFloeGroup::floe_type::point_type;

    PeriodicProblem() : base_class() {}
    PeriodicProblem(TSpaceTopology& topology) : base_class(), m_space_topology{topology} {}

    //! sets space borders adapted to floes limit positions
    void auto_topology();

private:
    TSpaceTopology m_space_topology;

    void set_topology_ptr(){
        base_class::m_dynamics_manager.set_topology(m_space_topology);
        base_class::m_proximity_detector.m_detector_h.set_topology(m_space_topology);
    }

};


template <
    typename TFloeGroup,
    typename TProxymityDetector,
    typename TCollisionManager,
    typename TDynamicsManager,
    typename TDomain,
    typename TSpaceTopology
>
void
PeriodicProblem<TFloeGroup, TProxymityDetector, TCollisionManager, TDynamicsManager, TDomain, TSpaceTopology>::
auto_topology()
{
    value_type min_x, min_y, max_x, max_y;
    min_x = min_y = std::numeric_limits<value_type>::max();
    max_x = max_y = - std::numeric_limits<value_type>::max();

    for (auto const& floe : base_class::m_floe_group.get_floes())  
        for (auto const& pt : floe.geometry().outer())
        {
            min_x = std::min(min_x, pt.x);
            min_y = std::min(min_y, pt.y);
            max_x = std::max(max_x, pt.x);
            max_y = std::max(max_y, pt.y);
        }

    value_type margin = 1;
    m_space_topology = TSpaceTopology{min_x - margin, max_x + margin, min_y - margin, max_y + margin};
    set_topology_ptr();
}



}} // namespace floe::problem


#endif // PROBLEM_PERIODIC_PROBLEM_HPP