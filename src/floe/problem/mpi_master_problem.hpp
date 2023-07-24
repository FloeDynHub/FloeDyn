/*!
 * \file floe/problem/mpi_master_problem.hpp
 * \brief Smooth mpi problem
 * \author Quentin Jouet
 */

#ifndef PROBLEM_MPI_MASTER_PROBLEM_HPP
#define PROBLEM_MPI_MASTER_PROBLEM_HPP

#include <iostream> // debug
#include <set>
#include <chrono> // tests

#include "floe/problem/mpi_problem.hpp"

namespace floe { namespace problem
{

/*! Problem
 *
 * It represents the whole problem of moving N floes in interval time [0, T], for concurrent runtime with mpi.
 *
 * \tparam TFloeGroup  
 * \tparam TProxymityDetector 
 * \tparam TCollisionManager   
 * \tparam TDynamicsManager 
 * \tparam TDomain
 *
 */

template <
    typename TProblem
>
class MPIMasterProblem : public MPIProblem<TProblem>
{
public:
    using base_class = MPIProblem<TProblem>;
    using real_type = typename TProblem::floe_group_type::floe_type::real_type;
    using point_type = typename TProblem::floe_group_type::floe_type::point_type;
    using message_type = typename base_class::message_type;
    using floe_distrib_type = typename TProblem::proximity_detector_type::floe_distrib_type;
    using process_list_type = typename TProblem::proximity_detector_type::process_list_type;

    //! Default constructor
    MPIMasterProblem(real_type epsilon, int OBL_status) : base_class(epsilon, OBL_status), msg_pk{0} {
        // this->m_out_manager.restrain_floe_ids({1, 2, 3});
    }

    //! Solver of the problem (main method)
    virtual void solve(real_type end_time, real_type dt_default, real_type out_step = 0, bool reset = true) override;
    virtual void recover_states_from_file(std::string const& filename, real_type t, bool keep_as_outfile=true) override;

private:
    //! last message id (increment for unicity)
    int msg_pk = 0; // TOD
    //! Move one time step forward
    virtual void step_solve() override;
     //! Collision solving
    virtual int manage_collisions() override;
    //! Compute next time step
    virtual void compute_time_step() override;
    //! Apply smooth dynamics to floes and verify interpenetration
    virtual void safe_move_floe_group() override;

    std::set<int> request_jobs(floe::io::JobTag, process_list_type const&, bool interpene=false);
    bool handle_responses(std::set<int>&, floe::io::JobTag);

    void send_request(message_type&, int);
    message_type receive_response();

    virtual void test_perf();
};


template<typename TProblem>
void MPIMasterProblem<TProblem>::solve(real_type end_time, real_type dt_default, real_type out_step, bool reset) {
    if (reset) this->create_optim_vars();
    this->m_domain.set_default_time_step(dt_default);
    this->m_out_manager.set_out_step(out_step, this->m_domain.time());
    this->output_datas(); // Initial state out
    while (this->m_domain.time() < end_time)
    {   
        auto t_start = std::chrono::high_resolution_clock::now();
        step_solve();
        auto t_end = std::chrono::high_resolution_clock::now();
        std::cout << "Chrono Time STEP : " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << std::endl;
        if (*this->QUIT) break; // exit normally after SIGINT
    }
    std::cout << " NB STEPS : " << this->m_step_nb << std::endl;
    // request_jobs(floe::io::termination_signal, this->m_proximity_detector.floe_process_distribution());
    // request_jobs(floe::io::termination_signal, this->m_proximity_detector.border_floe_process_distribution());
    // request_jobs(floe::io::termination_signal, this->m_proximity_detector.crossing_floe_process_distribution());
    request_jobs(floe::io::termination_signal, this->m_proximity_detector.all_worker_processes());
}

template<typename TProblem>
void MPIMasterProblem<TProblem>::step_solve(){
    auto t_start = std::chrono::high_resolution_clock::now();
    this->m_proximity_detector.distribute_floes();
    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << "Chrono distribute_floes : " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << std::endl;
    
    // test_perf();

    t_start = std::chrono::high_resolution_clock::now();
    manage_collisions();
    t_end = std::chrono::high_resolution_clock::now();
    std::cout << "Chrono Collision job : " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << std::endl;

    t_start = std::chrono::high_resolution_clock::now();
    compute_time_step();
    t_end = std::chrono::high_resolution_clock::now();
    std::cout << "Chrono Time Step job : " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << std::endl;

    t_start = std::chrono::high_resolution_clock::now();
    safe_move_floe_group();
    t_end = std::chrono::high_resolution_clock::now();
    std::cout << "Chrono Move job : " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << std::endl;
    if (this->m_dynamics_manager.get_external_forces().get_physical_data().get_air_mode()==5 || this->m_dynamics_manager.get_external_forces().get_physical_data().get_air_mode()==6) { //!< only if the external forces is a vortex
        for (size_t i=0; i<this->m_dynamics_manager.get_external_forces().get_physical_data().get_nb_vortex(); ++i) {
            std::cout << "the vortex wind speed is: " << 
                this->m_dynamics_manager.get_external_forces().get_physical_data().get_vortex_wind_speed(i) 
                << std::endl;
        }
	}
    // Output
    this->output_datas();
    if (this->m_step_nb % 10 == 0) this->m_proximity_detector.display_floe_distrib();

    this->m_step_nb++;
}

template<typename TProblem>
void MPIMasterProblem<TProblem>::test_perf(){
    auto t_start = std::chrono::high_resolution_clock::now();
    // auto msg_ids = request_jobs(floe::io::test_job, this->m_proximity_detector.floe_process_distribution());
    // auto msg_ids_2 = request_jobs(floe::io::test_job, this->m_proximity_detector.border_floe_process_distribution());
    // auto msg_ids_3 = request_jobs(floe::io::test_job, this->m_proximity_detector.crossing_floe_process_distribution());
    auto msg_ids = request_jobs(floe::io::test_job, this->m_proximity_detector.all_worker_processes());
    // msg_ids.insert(msg_ids_2.cbegin(), msg_ids_2.cend());
    // msg_ids.insert(msg_ids_3.cbegin(), msg_ids_3.cend());
    handle_responses(msg_ids, floe::io::test_job);
    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << "Chrono Empty job : " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << std::endl;
}

// template<typename TProblem>
// int MPIMasterProblem<TProblem>::manage_collisions(){
//     // TODO something more explicit !!
//     bool still_collision{true};
//     bool first{true};
//     int loop_count{0};
//     while (still_collision && loop_count < 20) {
//         loop_count++;
//         std::cout << "LCP : " << std::flush;
//         auto msg_ids = request_jobs(floe::io::collision_job, this->m_proximity_detector.floe_process_distribution());
//         still_collision = handle_responses(msg_ids, floe::io::collision_job);
//         if (!(first || still_collision)) break;
//         bool still_border_collision{true};
//         bool border_first{true};
//         int subloop_count{0};
//         while (still_border_collision && subloop_count <= 10){
//             subloop_count++;
//             std::cout << "Border LCP : " << std::flush;
//             msg_ids = request_jobs(floe::io::collision_job, this->m_proximity_detector.border_floe_process_distribution());
//             still_border_collision = handle_responses(msg_ids, floe::io::collision_job);
//             if (border_first) still_collision = still_border_collision;
//             if (!(border_first || still_border_collision)) break;
//             std::cout << "Cross LCP : " << std::flush;
//             msg_ids = request_jobs(floe::io::collision_job, this->m_proximity_detector.crossing_floe_process_distribution());
//             still_border_collision = handle_responses(msg_ids, floe::io::collision_job);
//             if (border_first) still_collision = still_collision || still_border_collision;
//             border_first = false;
//         }    
//         first = false;
//     } 
//     return 0;
// }

template<typename TProblem>
int MPIMasterProblem<TProblem>::manage_collisions(){
    // TODO something more explicit !!
    unsigned long OK_SET = 0; // first loop needs to consult all workers
    bool still_collision;
    unsigned long loop_count{0};
    auto keys = this->m_proximity_detector.collision_process_partition_keys();
    int i = 0;
    while (OK_SET < keys.size() and loop_count++ < 20 * keys.size()) {
        std::cout << "LCP (" << keys[i] << ") : " << std::flush;
        auto msg_ids = request_jobs(floe::io::collision_job, this->m_proximity_detector.process_partition().at(keys[i]));
        still_collision = handle_responses(msg_ids, floe::io::collision_job);
        if (!still_collision) { OK_SET++; } else { OK_SET = 0; }
        i = (i+1)%keys.size();
    }
    return 0;
}

template<typename TProblem>
void MPIMasterProblem<TProblem>::compute_time_step(){
    // auto msg_ids = request_jobs(floe::io::time_step_job, this->m_proximity_detector.floe_process_distribution());
    // auto msg_ids_2 = request_jobs(floe::io::time_step_job, this->m_proximity_detector.border_floe_process_distribution());
    // auto msg_ids_3 = request_jobs(floe::io::time_step_job, this->m_proximity_detector.crossing_floe_process_distribution());
    // msg_ids.insert(msg_ids_2.cbegin(), msg_ids_2.cend());
    // msg_ids.insert(msg_ids_3.cbegin(), msg_ids_3.cend());
    auto msg_ids = request_jobs(floe::io::time_step_job, this->m_proximity_detector.all_worker_processes());
    real_type delta_t = std::numeric_limits<real_type>::max();
    while (msg_ids.size()){
        auto resp = receive_response();
        msg_ids.erase(resp.id());
        delta_t = std::min(resp.time_step(), delta_t);
        // if (resp.time_step() < 5) std::cout << "#" << resp.mpi_source() << " " << resp.time_step() << std::endl;
    }
    this->m_domain.set_time_step(delta_t);
    // std::cout << "MC : " << this->m_floe_group.mass_center() << std::endl;
}

template<typename TProblem>
void MPIMasterProblem<TProblem>::safe_move_floe_group(){
    bool interpene{false};
    this->m_floe_group.backup_step_states();
    do {
        if (interpene){
            std::cout << "MASTER INTER" << std::endl;
            this->m_floe_group.recover_previous_step_states();
            this->m_domain.rewind_time();
            this->m_domain.set_time_step(this->m_domain.time_step() / 5);
            if (this->m_domain.time_step() < this->m_domain.default_time_step() / 1e8)
            {   
                // Hack to bypass repeating interpenetrations...
                this->m_out_manager.flush();
                this->recover_states_from_file(this->m_out_manager.out_file_name(), this->m_domain.time() + 1);
                std::cout << "dt too small -> RECOVER STATES FROM OUT FILE" << std::endl;
                return;
            }
        }
        auto msg_ids = request_jobs(floe::io::move_job, this->m_proximity_detector.base_grid_processes(), interpene);
        interpene = handle_responses(msg_ids, floe::io::move_job);
        if (!interpene){
            // auto msg_ids_2 = request_jobs(floe::io::interpene_job, this->m_proximity_detector.process_partition());
            // auto msg_ids_3 = request_jobs(floe::io::interpene_job, this->m_proximity_detector.crossing_floe_process_distribution());
            // msg_ids_2.insert(msg_ids_3.cbegin(), msg_ids_3.cend());
            auto msg_ids_2 = request_jobs(floe::io::interpene_job, this->m_proximity_detector.borders_processes());
            // border workers don't move any floe, they only check for interpenetration
            interpene = handle_responses(msg_ids_2, floe::io::interpene_job);
        }
        this->m_domain.update_time();
    } while (interpene);
}

// template<typename TProblem>
// std::set<int> MPIMasterProblem<TProblem>::request_jobs(floe::io::JobTag tag, floe_distrib_type const& floe_distrib, bool interpene){
//     std::set<int> msg_id_set;
//     for (auto const& iter : floe_distrib)
//     {
//         message_type request{++msg_pk};
//         msg_id_set.insert(msg_pk);
//         request.set_tag(tag);
//         request.store_time(this->m_domain.time());
//         request.set_floe_ids(iter.second);
//         request.store_states_light(this->get_floe_group(), request.floe_ids(), iter.first);
//         request.interpenetration(interpene);
//         if (tag==floe::io::move_job) {
//             request.store_time_step(this->m_domain.time_step());
//         }
//         send_request(request, iter.first); // iter.first = process id
//     }
//     return msg_id_set;
// }

template<typename TProblem>
std::set<int> MPIMasterProblem<TProblem>::request_jobs(floe::io::JobTag tag, process_list_type const& process_ids, bool interpene){
    std::set<int> msg_id_set;
    for (int p_id : process_ids)
    {
        message_type request{++msg_pk};
        msg_id_set.insert(msg_pk);
        request.set_tag(tag);
        request.store_time(this->m_domain.time());
        request.store_OBL_contribution(this->get_dynamics_manager().OBL_speed());
        request.set_floe_ids(this->m_proximity_detector.floe_process_distribution().at(p_id));
        request.store_states_light(this->get_floe_group(), request.floe_ids(), p_id);
        request.interpenetration(interpene);
        if (tag==floe::io::move_job) {
            request.store_time_step(this->m_domain.time_step());
        }
        send_request(request, p_id);
    }
    return msg_id_set;
}

template<typename TProblem>
bool MPIMasterProblem<TProblem>::handle_responses(std::set<int>& msg_id_set, floe::io::JobTag tag){
    bool ret{false};
    int lcp_tot = 0;
    point_type OBL_floes_force{0,0};
    while (msg_id_set.size()){ // while there is still workers who did not respond yet! 
        auto resp = receive_response();
        msg_id_set.erase(resp.id());
        this->get_floe_group().update_floe_states(resp, false); // do not update floe border or mesh (slow and useless for master)
        if (tag==floe::io::collision_job){ // TODO maybe use bool refs instead
            ret = ret || resp.nb_LCP_solved();
            int nb_lcp = resp.nb_LCP_solved();
            lcp_tot += nb_lcp;
            if (nb_lcp) std::cout << nb_lcp << ((msg_id_set.size()!=0) ? " + " : "") << std::flush;
            if (msg_id_set.size() == 0) std::cout << " = " << lcp_tot << std::endl;
        } else if (tag==floe::io::move_job or tag==floe::io::interpene_job) {
            OBL_floes_force += resp.template get_OBL_speed<point_type>();
            ret = ret || resp.interpenetration();
        }
    }
    if (tag==floe::io::move_job){
        this->get_dynamics_manager().update_ocean(this->get_floe_group(), this->m_domain.time_step(), OBL_floes_force);
    }
    return ret;
}

template<typename TProblem>
void MPIMasterProblem<TProblem>::send_request(message_type& request, int process_id){
    this->mpi().send_serial(request, process_id, 0);
}

template<typename TProblem>
typename MPIMasterProblem<TProblem>::message_type MPIMasterProblem<TProblem>::receive_response(){
    return this->mpi().template receive_serial<message_type>(MPI_ANY_SOURCE, MPI_ANY_TAG);
}

template<typename TProblem>
void MPIMasterProblem<TProblem>::recover_states_from_file(std::string const& filename, real_type t, bool keep_as_outfile){
    base_class::recover_states_from_file(filename, t);
    this->get_floe_group().post_load_floe(); // erase erroneous states origin
}

}} // namespace floe::problem


#endif // PROBLEM_MPI_MASTER_PROBLEM_HPP
