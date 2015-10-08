/*!
 * \file problem/periodic_problem.hpp
 * \brief Smooth periodic problem
 * \author Quentin Jouet
 */

#ifndef PROBLEM_PERIODIC_PROBLEM_HPP
#define PROBLEM_PERIODIC_PROBLEM_HPP

#include "floe/problem/problem.hpp"
#include "floe/io/matlab/pze_import.hpp"

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

    virtual inline void load_matlab_config(std::string const& filename) override {
        base_class::load_matlab_config(filename);
        auto imported_topo = floe::io::matlab::topology_from_file<TSpaceTopology>(filename);
        if (imported_topo.area() != 0)
            set_topology(imported_topo);
        else
            auto_topology();
    }

    //! sets space borders adapted to floes limit positions
    void set_topology(TSpaceTopology const& topology);
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
set_topology(TSpaceTopology const& topology)
{
    m_space_topology = topology;
    set_topology_ptr();
}


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
    auto a = base_class::m_floe_group.bounding_window();
    set_topology( TSpaceTopology{ a[0], a[1], a[2], a[3]} );
}



}} // namespace floe::problem


#endif // PROBLEM_PERIODIC_PROBLEM_HPP