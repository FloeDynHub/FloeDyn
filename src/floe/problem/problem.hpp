/*!
 * \file problem/problem.hpp
 * \brief Smooth problem
 * \author Quentin Jouet
 */

#ifndef PROBLEM_PROBLEM_HPP
#define PROBLEM_PROBLEM_HPP

#include "floe/problem/problem_h.hpp"
#include "floe/problem/problem_alg.hpp"
// #include "floe/io/false_hdf5_manager.hpp" // for gcc/MacOS
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
    using proximity_detector_type = TProxymityDetector;
    using floe_group_type = TFloeGroup;
    using collision_manager_type = TCollisionManager;
    using dynamics_manager_type = TDynamicsManager;
    using domain_type = TDomain;
    using problem_h_type = floe::problem::Problem_h<
        typename floe::problem::Problem_alg<
            typename floe_group_type::floe_group_h_type::floe_group_alg_type,
            typename collision_manager_type::manager_h_type::solver_type
        >,
        domain_type,
        typename floe_group_type::floe_group_h_type,
        typename proximity_detector_type::detector_h_type,
        typename collision_manager_type::manager_h_type
    >;
    using out_manager_type = io::HDF5Manager<floe_group_type, dynamics_manager_type>;

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

    virtual inline void load_matlab_config(std::string const& filename) {
        m_floe_group.load_matlab_config(filename);
        m_dynamics_manager.load_matlab_ocean_window_data(filename);
    }

    inline void load_matlab_topaz_data(std::string const& filename) {
        m_dynamics_manager.load_matlab_topaz_data(filename);
    }

    void recover_states_from_file(std::string const& filename, double t){
        double saved_time = m_out_manager.recover_states(filename, t, m_floe_group, m_dynamics_manager);
        m_domain.set_time(saved_time);
    }

    //! Solver
    void solve(double end_time, double out_step = 0){ // TODO change double
        create_optim_vars();
        m_domain.set_out_step(out_step);
        while (m_domain.time() < end_time)
        {
            step_solve();
            if (QUIT) break; // exit normally after SIGINT
        }
        std::cout << " NB STEPS : " << m_step_nb << std::endl;
    }

    void test(); // implemented in test file
    inline proximity_detector_type const& get_proximity_detector() const { return m_proximity_detector; }


protected:

    problem_h_type m_problem_h;

    // domain
    domain_type m_domain;

    // operators
    proximity_detector_type m_proximity_detector;
    collision_manager_type m_collision_manager;
    dynamics_manager_type m_dynamics_manager;

    // variables
    floe_group_type m_floe_group;

    int m_step_nb;
    out_manager_type m_out_manager;

    void create_optim_vars() {
        // mixes smooth and discrete levels because of detector structure
        // TODO improve this
        for (auto& floe_ptr : m_floe_group.get_floes())
            m_proximity_detector.m_detector_h.push_back(&floe_ptr);
    }

    void create_h(){
        m_problem_h.set_floe_group_h(m_floe_group.get_floe_group_h());
        m_problem_h.set_detector_h(m_proximity_detector.m_detector_h);
        m_problem_h.set_collision_manager_h(m_collision_manager.get_manager_h());
        m_problem_h.set_domain_h(m_domain);
        m_problem_h.create_alg();
    }

    // move one time step forward
    void step_solve(){
        m_problem_h.solve();
        m_dynamics_manager.move_floes(m_floe_group, m_domain.time_step());

        m_domain.update_time();
        std::cout << " Time : " << m_domain.time(); // << std::endl;
        std::cout << " | delta_t : " << m_domain.time_step();
        std::cout << " | Kinetic energy : " << m_floe_group.kinetic_energy() << std::endl;
        std::cout << "----" << std::endl;

        // ouput data
        if (m_domain.need_step_output())
        {
            m_out_manager.save_step(m_domain.time(), m_floe_group, m_dynamics_manager);
            m_domain.update_last_out();
            m_domain.update_next_out_limit();
        }

        m_step_nb++;
    }


};

}} // namespace floe::problem


#endif // PROBLEM_PROBLEM_HPP