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
 * \tparam TDomain_h            Discrete domain type
 * \tparam TFloeGroup_h         Discrete floe group variable type
 * \tparam TDetector            Detector type for collisions
 * \tparam TCollisionManager_h  Discrete collision manager type
 *
 */
template <
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

    using detector_h_type = TDetector;
    using time_scale_manager_type = ope::TimeScaleManager<detector_h_type>;


    //! Solver
    void solve() {
         step_solve();
    };

    //! Discrete floes group accessor
    inline TFloeGroup_h& get_floe_group_h() { return *m_floe_group_h; }
    inline void set_floe_group_h( TFloeGroup_h& floe_group_h){ m_floe_group_h = &floe_group_h; }

    inline void set_detector_h(detector_h_type& detector){ m_detector = &detector; }

    inline void set_collision_manager_h ( TCollisionManager_h& manager_h ) { 
        m_collision_manager_h = &manager_h; 
    }
    inline void set_domain_h(TDomain_h& domain){ m_domain_h = &domain; }

private:

    // domain
    TDomain_h* m_domain_h; //!< Discrete domain (same as smooth domain at the moment)

    // variable
    TFloeGroup_h* m_floe_group_h; //!< Discrete floes group (group of discrete floes)

    // operators
    detector_h_type* m_detector; //!< Proximity Detector at discrete level
    TCollisionManager_h* m_collision_manager_h; //!< Collision manager at discrete level
    time_scale_manager_type m_time_scale_manager; //!< Time scale manager at discrete level

    //! Proximity detection (inter-floe distance and eventual collisions)
    void do_detection(){
        if (!m_detector->update()) // we have a floe interpenetration
            m_domain_h->rewind_time();
    }

    //! Collision solving
    void manage_collisions(){
        m_collision_manager_h->solve_contacts(m_detector->contact_graph());
    }

    //! move one time step forward
    void step_solve(){
        do_detection();
        manage_collisions();
        m_time_scale_manager.delta_t_secu(m_domain_h, m_detector);
    }

    /* // chrono version (dev)
    void step_solve(){
        using namespace std;
        auto t_start = chrono::high_resolution_clock::now();
        do_detection();
        auto t_end = chrono::high_resolution_clock::now();
        cout << "detection : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;

        t_start = chrono::high_resolution_clock::now();
        manage_collisions();
        t_end = chrono::high_resolution_clock::now();
        cout << "collisions : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;

        t_start = chrono::high_resolution_clock::now();
        m_domain_h->set_time_step(m_time_scale_manager.delta_t_secu(
            m_domain_h, m_detector));
        t_end = chrono::high_resolution_clock::now();
        cout << "delta_t_secu : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;
    }*/


};

}} // namespace floe::problem


#endif // PROBLEM_PROBLEM_H_HPP