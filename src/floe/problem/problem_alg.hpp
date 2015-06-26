/*!
 * \file problem/problem.hpp
 * \brief Algebraic problem
 * \author Quentin Jouet
 */

#ifndef PROBLEM_PROBLEM_ALG_HPP
#define PROBLEM_PROBLEM_ALG_HPP


#include "floe/floes/static_floe.hpp"
#include "floe/floes/floe_exception.hpp"


namespace floe { namespace problem
{

/*! Problem_alg
 *
 * It represents the Step Algebraic problem.
 *
 * \tparam TFloeConfig_alg       Algebraic floe configuration type
 * \tparam TDomain_h             Discrete domain type
 * \tparam TExternalForces_alg   Algebraic external forces type
 * \tparam TCollisionSolver      Collision solver type
 * \tparam TForceIntegrator_alg  Algebraic force integrator type
 *
 */
template <
    typename TFloeConfig_alg,
    typename TDomain_h,
    typename TExternalForces_alg,
    typename TCollisionSolver,
    typename TForceIntegrator_alg,
>
class Problem_alg
{

public:

    //! Default constructor.
    // Problem_alg()

    //! Solver
    auto solve() {
        // LCP solve
        for (auto lcp : m_floe_config_alg->m_list_LCP){
            auto solution = m_collision_solver.solve_lcp(lcp);
            // update floes (maybe LCPSolver can do it)

        }
        // Fext integration
        m_force_integrator_alg->integrate(m_floe_config_alg->list_floe_alg,
                                          m_external_forces_alg,
                                          m_domain->m_time_delta)
        // update floes (a priori integrator does it)
    };

private:

    // domain
    TDomain_h* m_domain;

    // variable
    TFloeConfig_alg* m_floe_config_alg;
    TExternalForces_alg* m_external_forces_alg;

    // operators
    TCollisionSolver* m_collision_solver;
    TForceIntegrator_alg* m_force_integrator_alg;

};

}} // namespace floe::problem


#endif // PROBLEM_PROBLEM_ALG_HPP