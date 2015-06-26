/*!
 * \file problem/problem.hpp
 * \brief Discrete problem
 * \author Quentin Jouet
 */

#ifndef PROBLEM_PROBLEM_H_HPP
#define PROBLEM_PROBLEM_H_HPP


#include "floe/floes/static_floe.hpp"
#include "floe/floes/floe_exception.hpp"


namespace floe { namespace problem
{

/*! Problem_h
 *
 * It represents the discret problem.
 *
 * \tparam TProblem_alg         Associated algebraic problem type
 * \tparam TDomain_h            Discrete domain type
 * \tparam TFloeConfig_h        Discrete floe config variable type
 * \tparam TExternalForces_h    Discrete external Forces type
 * \tparam TDetector            Detector type for collisions
 * \tparam TCollisionManager_h  Discrete collision manager type
 * \tparam TForceIntegrator_h   Discrete force integrator type
 *
 */
template <
    typename TProblem_alg,
    typename TDomain_h,
    typename TFloeConfig_h,
    typename TExternalForces_h,
    typename TDetector,
    typename TCollisionManager_h,
    typename TForceIntegrator_h,
>
class Problem_h
{

public:

    //! Default constructor.
    Problem_h() : m_domain_h{nullptr}, m_external_forces_h{nullptr},
                  m_detector{nullptr}, m_collision_manager_h{nullptr},
                  m_force_integrator_h{nullptr} {}

    //! Solver
    auto solve() {
        while (domain_h.current_time < domain_h.domain.end_time){
            step_solve(domain_h.delta_t);
        };
    };

    // move one time step forward
    auto step_solve(){
        //detector
        //LCP creation
        m_problem_alg.solve()
    };

private:

    TProblem_alg m_problem_alg;

    // domain
    TDomain_h* m_domain_h;

    // variable
    TFloeConfig_h* m_floe_config_h;
    TExternalForces_h* m_external_forces_h;

    // operators
    TDetector* m_detector;
    TCollisionManager_h* m_collision_manager_h;
    TForceIntegrator_h* m_force_integrator_h;



};

}} // namespace floe::problem


#endif // PROBLEM_PROBLEM_H_HPP