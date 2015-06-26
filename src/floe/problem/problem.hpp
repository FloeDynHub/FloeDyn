/*!
 * \file problem/problem.hpp
 * \brief Smooth problem
 * \author Quentin Jouet
 */

#ifndef PROBLEM_PROBLEM_HPP
#define PROBLEM_PROBLEM_HPP

#include "floe/variable/floe_group.hpp"
#include "floe/ope/proximity_detector.hpp"
#include "floe/collision/matlab/detector.hpp"
// #include "floe/problem/problem_h.hpp"



//TODO
// #include "floe/ope/proxy-detector "
//TODO

namespace floe { namespace problem
{

/*! Problem
 *
 * It represents the whole problem of moving N floes in interval time [0, T].
 *
 * \tparam TFloe  
 * \tparam TDetector 
 * \tparam TProxymityDetector   
 *
 */

template <
    typename TFloe,
    typename TProxymityDetector
>
class Problem
{

public:

    // Type traits
    typedef TProxymityDetector          proximity_detector_type;
    // typedef TProblem_h                  problem_h_type;
    using floe_group_type = floe::variable::FloeGroup<TFloe>;

    //! Default constructor.
    // Problem()

    inline void load_matlab_config(std::string filename) {
        m_floe_group.load_matlab_config(filename);
    };

    inline void solve(){
        create_optim_vars();
    }

    //! Solver
    // auto solve() {
    //     m_problem_h.solve();
    // };

private:

    // problem_h_type m_problem_h;

    // domain
    // domain_type m_domain;

    // operators
    proximity_detector_type m_proximity_detector;
    // collision_manager_type m_collision_manager;
    // force_integrator_type m_force_integrator;

    // variables
    floe_group_type m_floe_group;
    // external_forces_type m_external_forces;

    inline void create_optim_vars() {
        for (auto& floe_ptr : m_floe_group.m_list_floe)
            m_proximity_detector.m_detector_h.push_back(&floe_ptr);
    };


};

}} // namespace floe::problem


#endif // PROBLEM_PROBLEM_HPP