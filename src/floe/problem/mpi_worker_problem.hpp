/*!
 * \file floe/problem/mpi_worker_problem.hpp
 * \brief Smooth mpi problem
 * \author Quentin Jouet
 */

#ifndef PROBLEM_MPI_WORKER_PROBLEM_HPP
#define PROBLEM_MPI_WORKER_PROBLEM_HPP

#include "floe/io/matlab/pze_import.hpp"

#include <iostream> // debug

#include "floe/problem/mpi_problem.hpp"
#include "floe/collision/contact_graph.hpp" // for graph vertices access (collision job response), todo move elsewhere

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
class MPIWorkerProblem : public MPIProblem<TProblem>
{
public:
    using base_class = MPIProblem<TProblem>;
    using real_type = typename TProblem::floe_group_type::floe_type::real_type;
    using point_type = typename TProblem::floe_group_type::floe_type::point_type;
    using message_type = typename base_class::message_type;

    //! Default constructor
    MPIWorkerProblem(real_type epsilon, int OBL_status) : base_class(epsilon, OBL_status), m_terminate{false} {}

    //! Solver of the problem (main method)
    virtual void solve(real_type end_time, real_type dt_default, real_type out_step = 0, bool reset = true, bool fracture = false, bool melting = false) override;

private:
    bool m_terminate;
    //! Move one time step forward
    virtual void step_solve(bool crack = false) override;
    message_type receive_request();
    void send_response(message_type&, message_type&);
};


template<typename TProblem>
void MPIWorkerProblem<TProblem>::solve(real_type end_time, real_type dt_default, real_type out_step, bool reset, bool fracture, bool melting) {
    if (reset) this->create_optim_vars();
    while (!m_terminate)
    {
        step_solve();
    }
}

template<typename TProblem>
void MPIWorkerProblem<TProblem>::step_solve(bool crack){
    // auto t_start = std::chrono::high_resolution_clock::now();
    message_type request = receive_request();
    // auto t_1 = std::chrono::high_resolution_clock::now();
    this->m_domain.set_time(request.time());
    this->get_floe_group().update_partial_list(request.floe_ids());
    this->get_floe_group().update_floe_states(request, true);
    for (auto& f : this->get_floe_group().get_floes()){
        // MPI + PBC has interpenetrations without this update (not sure exactly why as update_floe_states should do the job)
        f.update();
    }
    this->update_optim_vars(); // Needed for MPI + Periodic
    message_type response{request};
    // auto t_10 = std::chrono::high_resolution_clock::now();
    if (request.tag()==floe::io::collision_job){
        if (!this->m_proximity_detector.update()) std::cout << "DIRECT INTER #" << this->mpi().process_rank() << std::endl;
        int n = this->manage_collisions();
        response.nb_LCP_solved(n);
    } else if (request.tag()==floe::io::time_step_job) {
        this->compute_time_step();
        response.store_time_step(this->m_domain.time_step());
    } else if (request.tag()==floe::io::move_job) {
        this->m_domain.set_time_step(request.time_step());
        this->get_dynamics_manager().set_OBL_speed(request.template get_OBL_speed<point_type>());
        point_type diff_speed = this->move_floe_group();
        bool interpene = !this->m_proximity_detector.check_interpenetration();
        response.interpenetration(interpene);
        response.store_OBL_contribution(diff_speed);
    } else if (request.tag()==floe::io::interpene_job){
        response.interpenetration(!this->m_proximity_detector.check_interpenetration());
    } else if (request.tag()==floe::io::termination_signal){
        m_terminate = true;
        return;
    }
    this->m_step_nb++;
    // auto t_2 = std::chrono::high_resolution_clock::now();
    send_response(response, request);
    // auto t_end = std::chrono::high_resolution_clock::now();
    // std::cout << "Chrono WORKER (job " << request.tag() << " on #" << this->mpi().process_rank() << ") "
        // << std::chrono::duration<double, std::milli>(t_10-t_1).count() << " + " 
        // << std::chrono::duration<double, std::milli>(t_2-t_10).count() << " + "
        // << std::chrono::duration<double, std::milli>(t_end-t_2).count() << " = "
        // << std::chrono::duration<double, std::milli>(t_end-t_1).count() << " ms"
        // << std::endl;
}

template<typename TProblem>
typename MPIWorkerProblem<TProblem>::message_type MPIWorkerProblem<TProblem>::receive_request(){
    return this->mpi().template receive_serial<message_type>(0, MPI_ANY_TAG);
}

template<typename TProblem>
void MPIWorkerProblem<TProblem>::send_response(message_type& resp, message_type& req){
    if (resp.tag() != floe::io::time_step_job){
        if (resp.tag() == floe::io::collision_job){

        } else {
            resp.store_states(this->get_floe_group(), req.floe_ids());
        }
        
    }
    switch(resp.tag()) {
        case floe::io::time_step_job : break;
        case floe::io::interpene_job : break;
        case floe::io::collision_job : {
            // if no LCP were solved, no floe states were modified
            if (resp.nb_LCP_solved() > 0){
                // Send only states of potentially changed floes (floes in contact graph)
                // TODO -> Useless as is because all floes are vertex of the contact graph
                //      -> collect from edges ?
                auto const& contact_graph = this->m_proximity_detector.contact_graph();
                std::vector<std::size_t> ids;
                for ( auto const v : boost::make_iterator_range( vertices(contact_graph)) ){
                    ids.push_back(contact_graph[v].floe->id());
                }
                resp.store_states(this->get_floe_group(), ids);
            }
            break;
        }
        default : resp.store_states(this->get_floe_group(), req.floe_ids()); break;
    }


    this->mpi().send_serial(resp, 0, 0);
}


}} // namespace floe::problem


#endif // PROBLEM_MPI_WORKER_PROBLEM_HPP