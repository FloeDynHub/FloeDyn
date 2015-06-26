/*!
 * \file problem/problem.hpp
 * \brief Smooth problem
 * \author Quentin Jouet
 */

#ifndef PROBLEM_PROBLEM_HPP
#define PROBLEM_PROBLEM_HPP

#include "floe/problem/problem_h.hpp"
#include "floe/problem/problem_alg.hpp"

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

    //! Default constructor.
    Problem() : m_step_nb{0}
    {
        create_h();
    }

    inline void load_matlab_config(std::string filename) {
        m_floe_group.load_matlab_config(filename);
    }

    //! Solver
    inline void solve(double end_time, double out_step = 0){ // TODO change double
        create_optim_vars();
        while (m_domain.time() < end_time)
            step_solve(out_step);
    }

    void test(); // implemented in test file
    inline proximity_detector_type const& get_proximity_detector() const { return m_proximity_detector; }


private:

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

    inline void create_optim_vars() {
        // mixes smooth and discrete levels because of detector structure
        // TODO improve this
        for (auto& floe_ptr : m_floe_group.get_floes())
            m_proximity_detector.m_detector_h.push_back(&floe_ptr);
    }

    inline void create_h(){
        m_problem_h.set_floe_group_h(m_floe_group.get_floe_group_h());
        m_problem_h.set_detector_h(m_proximity_detector.m_detector_h);
        m_problem_h.set_collision_manager_h(m_collision_manager.get_manager_h());
        m_problem_h.set_domain_h(m_domain);
        m_problem_h.create_alg();
    }

    // move one time step forward
    inline void step_solve(double out_step = 0){
        m_problem_h.solve();
        m_dynamics_manager.move_floes(m_floe_group, m_domain.time_step());

        m_domain.update_time();
        std::cout << " Time : " << m_domain.time() << std::endl << "----" << std::endl;

        // if (out_step && m_step_nb % out_step == 0)
        //     m_floe_group.out_hdf5(m_domain.time());

        if (out_step && m_domain.time() - m_domain.last_out() >= out_step)
        {
            m_floe_group.out_hdf5(m_domain.time());
            m_domain.update_last_out();
        }

        m_step_nb++;
    }


};

}} // namespace floe::problem


#endif // PROBLEM_PROBLEM_HPP