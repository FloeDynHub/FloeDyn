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
    using real_type = typename TFloeGroup::floe_type::real_type;
    using point_type = typename TFloeGroup::floe_type::point_type;

    //! Default constructor
    PeriodicProblem(real_type epsilon, int OBL_status) : base_class(epsilon, OBL_status) {}
    //! Constructor from topology
    PeriodicProblem(TSpaceTopology& topology) : base_class(), m_space_topology{topology} {}

    //! Load initial state and construct topology from matlab file
    virtual void load_matlab_config(std::string const& filename) override;

    //! construct topology from box information in Matlab file
    void load_topology_from_matlab(std::string const& filename);

    //! Load initial state and construct topology from matlab file
    virtual void load_h5_config(std::string const& filename) override;

    //! Set topology
    void set_topology(TSpaceTopology const& topology);
    //! Sets topology calculating space borders adapted to floes limit positions
    void auto_topology();
    //! Floe Concentration override (takes topology window area as reference)
    virtual real_type floe_concentration() override { return base_class::m_floe_group.total_area() / m_space_topology.area(); }

private:
    TSpaceTopology m_space_topology;

    void set_topology_ptr(){
        base_class::m_dynamics_manager.set_topology(m_space_topology);
        if (base_class::m_proximity_detector.is_periodic()) {
            base_class::m_proximity_detector.set_topology(m_space_topology);
        }
    }

};

// MACRO def to lighten file
#define TEMPLATE_PERIO_PB template <\
    typename TFloeGroup,\
    typename TProxymityDetector,\
    typename TCollisionManager,\
    typename TDynamicsManager,\
    typename TDomain,\
    typename TSpaceTopology\
>
#define PERIODIC_PROBLEM PeriodicProblem<TFloeGroup, TProxymityDetector, TCollisionManager, TDynamicsManager, TDomain, TSpaceTopology>


TEMPLATE_PERIO_PB
void PERIODIC_PROBLEM::load_matlab_config(std::string const& filename) {
    base_class::load_matlab_config(filename);
    load_topology_from_matlab(filename);
}


TEMPLATE_PERIO_PB
void PERIODIC_PROBLEM::load_topology_from_matlab(std::string const& filename)  {
    auto imported_topo = floe::io::matlab::topology_from_file<TSpaceTopology>(filename);
    if (imported_topo.area() != 0)
        set_topology(imported_topo);
    else
        auto_topology();
}


TEMPLATE_PERIO_PB
void PERIODIC_PROBLEM::load_h5_config(std::string const& filename) {
    base_class::load_h5_config(filename);
    if (base_class::m_floe_group.initial_window_area())
    {
        auto a = base_class::m_floe_group.get_initial_window();
        set_topology( TSpaceTopology{ a[0], a[1], a[2], a[3]} );
    } else
        auto_topology();
}


TEMPLATE_PERIO_PB
void PERIODIC_PROBLEM::set_topology(TSpaceTopology const& topology)
{
    m_space_topology = topology;
    set_topology_ptr();
}


TEMPLATE_PERIO_PB
void PERIODIC_PROBLEM::auto_topology()
{
    auto a = base_class::m_floe_group.bounding_window();
    set_topology( TSpaceTopology{ a[0], a[1], a[2], a[3]} );
}



}} // namespace floe::problem


#endif // PROBLEM_PERIODIC_PROBLEM_HPP