/*!
 * \file problem/problem.hpp
 * \brief Smooth problem
 * \author Quentin Jouet
 */

#ifndef PROBLEM_PROBLEM_HPP
#define PROBLEM_PROBLEM_HPP

#include "floe/problem/problem_h.hpp"
#include "floe/io/hdf5_manager.hpp"

#include <iostream> // debug

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

    //! Default constructor.
    Problem() :
        m_problem_h{},
        m_domain{},
        m_proximity_detector{},
        m_collision_manager{},
        m_dynamics_manager{m_domain.time()},
        m_floe_group{},
        m_step_nb{0}
    {
        create_h();
    }

    //! Load floes set and initial states from matlab file
    virtual inline void load_matlab_config(std::string const& filename) {
        m_floe_group.load_matlab_config(filename);
        m_dynamics_manager.load_matlab_ocean_window_data(filename, m_floe_group);
    }

    //! Load ocean and wind data from a topaz file
    inline void load_matlab_topaz_data(std::string const& filename) {
        m_dynamics_manager.load_matlab_topaz_data(filename);
    }

    //! Recover simulation state from previous ouput file, at any recorded time t
    void recover_states_from_file(std::string const& filename, value_type t){
        value_type saved_time = m_out_manager.recover_states(filename, t, m_floe_group, m_dynamics_manager);
        m_domain.set_time(saved_time);
    }

    //! Solver of the problem (main method)
    void solve(value_type end_time, value_type dt_default, value_type OBL_status, value_type out_step = 0){
        create_optim_vars();
        m_domain.set_out_step(out_step);
        m_domain.set_default_time_step(dt_default);
        m_dynamics_manager.set_OBL_status(OBL_status);
        while (m_domain.time() < end_time)
        {
            step_solve();
            if (QUIT) break; // exit normally after SIGINT
        }
        std::cout << " NB STEPS : " << m_step_nb << std::endl;
    }


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
    void create_optim_vars() {
        // mixes smooth and discrete levels because of detector structure
        // TODO improve this
        for (auto& floe_ptr : m_floe_group.get_floes())
            m_proximity_detector.m_detector_h.push_back(&floe_ptr);
    }

    //! Initializing discrete problem
    void create_h(){
        m_problem_h.set_floe_group_h(m_floe_group.get_floe_group_h());
        m_problem_h.set_detector_h(m_proximity_detector.m_detector_h);
        m_problem_h.set_collision_manager_h(m_collision_manager.get_manager_h());
        m_problem_h.set_domain_h(m_domain);
    }

    //! Move one time step forward
    void step_solve(){
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
            m_out_manager.save_step(m_domain.time(), m_floe_group, m_dynamics_manager);
            // m_domain.update_last_out();
            m_domain.update_next_out_limit();
        }

        m_step_nb++;
    }


};

}} // namespace floe::problem


#endif // PROBLEM_PROBLEM_HPP