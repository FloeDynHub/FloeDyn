/*!
 * \file problem/problem.hpp
 * \brief Smooth problem
 * \author Quentin Jouet
 */

#ifndef PROBLEM_PROBLEM_HPP
#define PROBLEM_PROBLEM_HPP

#include "floe/problem/problem_h.hpp"
#include "floe/io/hdf5_manager.h"

#include <iostream> // debug
#include <atomic>

namespace floe { namespace problem
{

/*! Problem
 *
 * It represents the whole problem of moving N floes in interval time [0, T].
 *
 * \tparam TFloeGroup  
 * \tparam TProxymityDetector 
 * \tparam TCollisionManager   
 * \tparam TDynamicsManager 
 * \tparam TDomain
 *
 */

template <
    typename TFloeGroup,
    typename TProxymityDetector,
    typename TCollisionManager,
    typename TDynamicsManager,
    typename TDomain
>
class Problem
{

public:

    // Type traits
    using problem_h_type = floe::problem::Problem_h<
        TDomain,
        typename TFloeGroup::floe_group_h_type,
        typename TProxymityDetector::detector_h_type,
        typename TCollisionManager::manager_h_type
    >;
    using out_manager_type = io::HDF5Manager<TFloeGroup, TDynamicsManager>;
    using value_type = typename TFloeGroup::value_type;
    using floe_group_type = TFloeGroup; // generator accessor

    //! Default constructor.
    Problem();

    //! Load floes set and initial states from matlab file
    virtual void load_matlab_config(std::string const& filename);

    //! Load ocean and wind data from a topaz file
    inline void load_matlab_topaz_data(std::string const& filename) {
        m_dynamics_manager.load_matlab_topaz_data(filename);
    }

    //! Recover simulation state from previous ouput file, at any recorded time t
    void recover_states_from_file(std::string const& filename, value_type t);

    //! Solver of the problem (main method)
    void solve(value_type end_time, value_type dt_default, value_type OBL_status, value_type out_step = 0, bool reset = true);

    //! set existing floe_group (from generator for example)
    virtual void set_floe_group(floe_group_type& floe_group, value_type win_area = 0);
    //! Floe group accessor for config generator
    inline TFloeGroup& get_floe_group() { return m_floe_group; }
    //! Dynamics mgr accessor for config generator
    inline TDynamicsManager& get_dynamics_manager(){ return m_dynamics_manager; }

    //! Floe Concentration
    virtual value_type floe_concentration() { return m_floe_group.floe_concentration(); }

    const std::atomic<bool>* QUIT; //!< Exit signal

protected:

    problem_h_type m_problem_h; //!< Discrete problem

    // domain
    TDomain m_domain; //!< Domain (time manager at the moment)

    // operators
    TProxymityDetector m_proximity_detector; //!< Object managing floes proximity (distance and collision detection)
    TCollisionManager m_collision_manager; //!< Object managing collisions solving
    TDynamicsManager m_dynamics_manager; //!< Object managing floes dynamics (moving according to physics)

    // variables
    TFloeGroup m_floe_group; //!< The set of floes

    int m_step_nb; //!< Total number of steps from beginning
    out_manager_type m_out_manager; //!< Object managing simulation output

    //! Initializing proximity detector with floe set
    void create_optim_vars();

    //! Initializing discrete problem
    void create_h();

    //! Move one time step forward
    void step_solve();

};


template <
    typename TFloeGroup,
    typename TProxymityDetector,
    typename TCollisionManager,
    typename TDynamicsManager,
    typename TDomain
>
Problem<TFloeGroup, TProxymityDetector,TCollisionManager, TDynamicsManager, TDomain>
::Problem() :
        QUIT{nullptr},
        m_problem_h{},
        m_domain{},
        m_proximity_detector{},
        m_collision_manager{},
        m_dynamics_manager{m_domain.time()},
        m_floe_group{},
        m_step_nb{0},
        m_out_manager{m_floe_group}
    {
        create_h();
    }


template <
    typename TFloeGroup,
    typename TProxymityDetector,
    typename TCollisionManager,
    typename TDynamicsManager,
    typename TDomain
>
void Problem<TFloeGroup, TProxymityDetector,TCollisionManager, TDynamicsManager, TDomain>
::load_matlab_config(std::string const& filename) {
    m_floe_group.load_matlab_config(filename);
    m_dynamics_manager.load_matlab_ocean_window_data(filename, m_floe_group);
}


template <
    typename TFloeGroup,
    typename TProxymityDetector,
    typename TCollisionManager,
    typename TDynamicsManager,
    typename TDomain
>
void Problem<TFloeGroup, TProxymityDetector,TCollisionManager, TDynamicsManager, TDomain>
::recover_states_from_file(std::string const& filename, value_type t){
    value_type saved_time = m_out_manager.recover_states(filename, t, m_floe_group, m_dynamics_manager);
    m_domain.set_time(saved_time);
}


template <
    typename TFloeGroup,
    typename TProxymityDetector,
    typename TCollisionManager,
    typename TDynamicsManager,
    typename TDomain
>
void Problem<TFloeGroup, TProxymityDetector,TCollisionManager, TDynamicsManager, TDomain>
::solve(value_type end_time, value_type dt_default, value_type OBL_status, value_type out_step, bool reset){
    if (reset) create_optim_vars();
    m_domain.set_out_step(out_step);
    m_domain.set_default_time_step(dt_default);
    m_dynamics_manager.set_OBL_status(OBL_status);
    while (m_domain.time() < end_time)
    {
        step_solve();
        if (*QUIT) break; // exit normally after SIGINT
    }
    std::cout << " NB STEPS : " << m_step_nb << std::endl;
}


template <
    typename TFloeGroup,
    typename TProxymityDetector,
    typename TCollisionManager,
    typename TDynamicsManager,
    typename TDomain
>
void Problem<TFloeGroup, TProxymityDetector,TCollisionManager, TDynamicsManager, TDomain>
::set_floe_group(floe_group_type& floe_group, value_type win_area) {
    m_floe_group = std::move(floe_group);
    m_dynamics_manager.set_ocean_window_area((win_area != 0) ? win_area : m_floe_group.bounding_window_area());
    m_out_manager.set_floe_group(m_floe_group);
}


template <
    typename TFloeGroup,
    typename TProxymityDetector,
    typename TCollisionManager,
    typename TDynamicsManager,
    typename TDomain
>
void Problem<TFloeGroup, TProxymityDetector,TCollisionManager, TDynamicsManager, TDomain>
::create_optim_vars() {
    // mixes smooth and discrete levels because of detector structure
    // TODO improve this
    m_proximity_detector.m_detector_h.reset();
    for (auto& floe : m_floe_group.get_floes())
        m_proximity_detector.m_detector_h.push_back(&floe);
}


template <
    typename TFloeGroup,
    typename TProxymityDetector,
    typename TCollisionManager,
    typename TDynamicsManager,
    typename TDomain
>
void Problem<TFloeGroup, TProxymityDetector,TCollisionManager, TDynamicsManager, TDomain>
::create_h(){
    m_problem_h.set_floe_group_h(m_floe_group.get_floe_group_h());
    m_problem_h.set_detector_h(m_proximity_detector.m_detector_h);
    m_problem_h.set_collision_manager_h(m_collision_manager.get_manager_h());
    m_problem_h.set_domain_h(m_domain);
}


template <
    typename TFloeGroup,
    typename TProxymityDetector,
    typename TCollisionManager,
    typename TDynamicsManager,
    typename TDomain
>
void Problem<TFloeGroup, TProxymityDetector,TCollisionManager, TDynamicsManager, TDomain>
::step_solve(){
    m_problem_h.solve();
    m_dynamics_manager.move_floes(m_floe_group, m_domain.time_step());

    m_domain.update_time();
    std::cout << " Time : " << m_domain.time();
    std::cout << " | delta_t : " << m_domain.time_step();
    std::cout << " | Kinetic energy : " << m_floe_group.kinetic_energy() << std::endl;
    std::cout << "----" << std::endl;

    // ouput data
    if (m_domain.need_step_output())
    {
        m_out_manager.save_step(m_domain.time(), m_dynamics_manager);
        // m_domain.update_last_out();
        m_domain.update_next_out_limit();
    }

    m_step_nb++;
}


}} // namespace floe::problem


#endif // PROBLEM_PROBLEM_HPP