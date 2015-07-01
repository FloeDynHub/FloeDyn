/*!
 * \file problem/problem.hpp
 * \brief Discrete problem
 * \author Quentin Jouet
 */

#ifndef PROBLEM_PROBLEM_H_HPP
#define PROBLEM_PROBLEM_H_HPP


#include "floe/floes/static_floe.hpp"
#include "floe/floes/floe_exception.hpp"
#include "floe/lcp/builder/graph_to_lcp.hpp"
#include "floe/ope/time_scale_manager.hpp"

#include <iostream> // debug
#include <chrono>


namespace floe { namespace problem
{
// using namespace boost;

/*! Problem_h
 *
 * It represents the discrete problem.
 *
 * \tparam TProblem_alg         Associated algebraic problem type
 * \tparam TDomain_h            Discrete domain type
 * \tparam TFloeGroup_h         Discrete floe group variable type
 * \tparam TDetector            Detector type for collisions
 * \tparam TCollisionManager_h  Discrete collision manager type
 *
 */
template <
    typename TProblem_alg,
    typename TDomain_h,
    typename TFloeGroup_h,
    typename TDetector,
    typename TCollisionManager_h
>
class Problem_h
{

public:
     //! Default constructor.
    Problem_h() : m_floe_group_h{nullptr}, m_detector{nullptr}, 
                  m_collision_manager_h{nullptr} {}

    using real = double; // TODO do better
    using floe_group_h_type = TFloeGroup_h;
    using detector_h_type = TDetector;
    using time_scale_manager_type = ope::TimeScaleManager<TDomain_h, detector_h_type>;


    inline void create_alg(){
        m_problem_alg.set_floe_group_alg(m_floe_group_h->get_floe_group_alg());
        m_problem_alg.set_collision_solver(m_collision_manager_h->get_solver());
    }

    //! Solver
    void solve() {
         step_solve();
    };

    inline TFloeGroup_h& get_floe_group_h() { return *m_floe_group_h; }
    inline void set_floe_group_h( TFloeGroup_h& floe_group_h){ m_floe_group_h = &floe_group_h; }

    inline void set_detector_h(detector_h_type& detector){ m_detector = &detector; }

    inline void set_collision_manager_h ( TCollisionManager_h& manager_h ) { 
        m_collision_manager_h = &manager_h; 
    }
    inline void set_domain_h(TDomain_h& domain){ m_domain_h = &domain; }

private:

    TProblem_alg m_problem_alg;

    // domain
    TDomain_h* m_domain_h;

    // variable
    TFloeGroup_h* m_floe_group_h;

    // operators
    detector_h_type* m_detector;
    TCollisionManager_h* m_collision_manager_h;
    time_scale_manager_type m_time_scale_manager;

    void do_detection(){
        m_detector->update();
    }

    void manage_collisions(){
        m_collision_manager_h->solve_contacts(m_detector->contact_graph());
    }

    //! move one time step forward
    void step_solve(){
        do_detection();
        manage_collisions();
        m_domain_h->set_time_step(m_time_scale_manager.delta_t_secu(
            m_domain_h, m_detector));
    }

    // chrono version (dev)
    // void step_solve(){
    //     using namespace std;
    //     auto t_start = chrono::high_resolution_clock::now();
    //     do_detection();
    //     auto t_end = chrono::high_resolution_clock::now();
    //     cout << "detection : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;

    //     t_start = chrono::high_resolution_clock::now();
    //     manage_collisions();
    //     t_end = chrono::high_resolution_clock::now();
    //     cout << "collisions : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;

    //     t_start = chrono::high_resolution_clock::now();
    //     m_domain_h->set_time_step(m_time_scale_manager.delta_t_secu(
    //         m_domain_h, m_detector));
    //     t_end = chrono::high_resolution_clock::now();
    //     cout << "delta_t_secu : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;
    // }


};

}} // namespace floe::problem


#endif // PROBLEM_PROBLEM_H_HPP