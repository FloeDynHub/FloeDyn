/*!
 * \file problem/problem.hpp
 * \brief Algebraic problem
 * \author Quentin Jouet
 */

#ifndef PROBLEM_PROBLEM_ALG_HPP
#define PROBLEM_PROBLEM_ALG_HPP


#include "floe/floes/static_floe.hpp"
#include "floe/floes/floe_exception.hpp"
#include <iostream> // DEBUG


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
    typename TFloeGroup_alg,
    // typename TDomain_h,
    // typename TExternalForces_alg,
    typename TCollisionSolver
    // typename TForceIntegrator_alg,
>
class Problem_alg
{

public:

    using floe_group_alg_type = TFloeGroup_alg;

    //! Default constructor.
    Problem_alg() : m_floe_group_alg{nullptr}, m_collision_solver{nullptr} {}

    //! Solver
    bool solve() {
        // LCP solve
        std::cout << m_floe_group_alg->get_list_LCP().size() << " LCP :"; //DEBUG
        int nb_success = 0;
        for (auto& lcp : m_floe_group_alg->get_list_LCP()){
            auto success = m_collision_solver->solve(lcp);
            // std::cout << std::endl << "LCP SOLVE SUCCESS : " << success << std::endl; // DEBUG
            if (success){ nb_success++; }
            // update floes (maybe LCPSolver can do it)
        }

        // clearing LCP list
        // m_floe_group_alg->empty_list_LCP();

        std::cout << nb_success << " SUCCESS"; //DEBUG
        return (nb_success != 0);

    };

    inline TFloeGroup_alg& get_floe_group_alg(){ return *m_floe_group_alg; }
    inline void set_floe_group_alg( TFloeGroup_alg& floe_group_alg){ m_floe_group_alg = &floe_group_alg; }

    inline void set_collision_solver( TCollisionSolver& solver){ m_collision_solver = &solver; }

private:

    // domain
    // TDomain_h* m_domain;

    // variable
    TFloeGroup_alg* m_floe_group_alg;
    // TExternalForces_alg* m_external_forces_alg;

    // operators
    TCollisionSolver* m_collision_solver;
    // TForceIntegrator_alg* m_force_integrator_alg;

};

}} // namespace floe::problem


#endif // PROBLEM_PROBLEM_ALG_HPP