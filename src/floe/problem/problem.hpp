/*!
 * \file problem/problem.hpp
 * \brief Smooth problem
 * \author Quentin Jouet
 */

#ifndef PROBLEM_PROBLEM_HPP
#define PROBLEM_PROBLEM_HPP

#include "floe/io/hdf5_manager.h"
// #include "floe/io/false_hdf5_manager.hpp"

 #include "floe/domain/time_scale_manager.hpp"

#include <iostream>
#include <atomic>
#include <stdexcept>

 #include <numeric> // FOR TEST

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
    using out_manager_type = io::HDF5Manager<TFloeGroup, TDynamicsManager>;
    using value_type = typename TFloeGroup::value_type;
    using floe_group_type = TFloeGroup; // generator accessor
    using time_scale_manager_type = domain::TimeScaleManager<typename TProxymityDetector::proximity_data_type>;
    using proximity_detector_type = TProxymityDetector;

    //! Default constructor.
    Problem(value_type epsilon=0.4, int OBL_status=0);

    //! Solver of the problem (main method)
    virtual void solve(value_type end_time, value_type dt_default, value_type out_step = 0, bool reset = true);

    virtual void load_config(std::string const& filename);
    //! Load ocean and wind data from a topaz file
    inline void load_matlab_topaz_data(std::string const& filename) {
        m_dynamics_manager.load_matlab_topaz_data(filename);
    }
    //! Recover simulation state from previous ouput file, at any recorded time t
    virtual void recover_states_from_file(std::string const& filename, value_type t, bool keep_as_outfile=true);

    //! set existing floe_group (from generator for example)
    virtual void set_floe_group(floe_group_type& floe_group);
    //! Floe group accessor for config generator
    inline TFloeGroup& get_floe_group() { return m_floe_group; }
    //! Dynamics mgr accessor for config generator
    inline TDynamicsManager& get_dynamics_manager(){ return m_dynamics_manager; }

    //! Floe Concentration
    virtual value_type floe_concentration() { return m_floe_group.floe_concentration(); }

    void make_input_file();

    const std::atomic<bool>* QUIT; //!< Exit signal

protected:
    // domain
    TDomain m_domain; //!< Domain (time manager at the moment)

    // operators
    TProxymityDetector m_proximity_detector; //!< Object managing floes proximity (distance and collision detection)
    TCollisionManager m_collision_manager; //!< Object managing collisions solving
    TDynamicsManager m_dynamics_manager; //!< Object managing floes dynamics (moving according to physics)
    time_scale_manager_type m_time_scale_manager; //!< Time scale manager at discrete level

    // variables
    TFloeGroup m_floe_group; //!< The set of floes

    // io
    int m_step_nb; //!< Total number of steps from beginning
    out_manager_type m_out_manager; //!< Object managing simulation output

    //! Load floes set and initial states from matlab file
    virtual void load_matlab_config(std::string const& filename);
    //! Load floes set and initial states from hdf5 file
    virtual void load_h5_config(std::string const& filename);
    //! Initializing proximity detector with floe set
    void create_optim_vars();
    //! Move one time step forward
    virtual void step_solve();
    //! Proximity detection (inter-floe distance and eventual collisions)
    void detect_proximity();
    //! Collision solving
    virtual int manage_collisions();
    //! Compute next time step
    virtual void compute_time_step();
    //! Apply smooth dynamics to floes and verify interpenetration
    virtual void safe_move_floe_group();
    //! Apply smooth dynamics to floes
    void move_floe_group();
    //! Handle output_datas (console + out file)
    void output_datas();
};


// MACRO def to lighten file
#define TEMPLATE_PB template <\
    typename TFloeGroup,\
    typename TProxymityDetector,\
    typename TCollisionManager,\
    typename TDynamicsManager,\
    typename TDomain\
>
#define PROBLEM Problem<TFloeGroup, TProxymityDetector, TCollisionManager, TDynamicsManager, TDomain>


TEMPLATE_PB
PROBLEM::Problem(value_type epsilon, int OBL_status) :
        QUIT{nullptr},
        m_domain{},
        m_proximity_detector{},
        m_collision_manager{epsilon},
        m_dynamics_manager{m_domain.time(), OBL_status},
        m_floe_group{},
        m_step_nb{0},
        m_out_manager{m_floe_group}
    {
        m_time_scale_manager.set_prox_data_ptr( &(m_proximity_detector.data()) );
    }

TEMPLATE_PB
void PROBLEM::load_config(std::string const& filename) {
    std::string extension = filename.substr(filename.find(".") + 1);
    if (extension == "mat")
        load_matlab_config(filename);
    else if (extension == "h5")
        load_h5_config(filename);
    else
        throw std::runtime_error("Input file extension must be .mat or .h5");
}

TEMPLATE_PB
void PROBLEM::load_matlab_config(std::string const& filename) {
    m_floe_group.load_matlab_config(filename);
    m_dynamics_manager.load_matlab_ocean_window_data(filename, m_floe_group);
}

TEMPLATE_PB
void PROBLEM::load_h5_config(std::string const& filename) {
    m_floe_group.load_h5_config(filename);
    m_dynamics_manager.set_ocean_window_area(m_floe_group.ocean_window_area());
}


TEMPLATE_PB
void PROBLEM::recover_states_from_file(std::string const& filename, value_type t, bool keep_as_outfile){
    value_type saved_time = m_out_manager.recover_states(filename, t, m_floe_group, m_dynamics_manager, keep_as_outfile);
    m_domain.set_time(saved_time);
}


TEMPLATE_PB
void PROBLEM::set_floe_group(floe_group_type& floe_group) {
    m_floe_group = std::move(floe_group);
    m_dynamics_manager.set_ocean_window_area(m_floe_group.ocean_window_area());
    m_out_manager.set_floe_group(m_floe_group);
}


TEMPLATE_PB
void PROBLEM::create_optim_vars() {
    m_proximity_detector.reset();
    // for (auto& floe : m_floe_group.get_floes())
    //     m_proximity_detector.push_back(&floe);
    m_proximity_detector.set_floe_group(m_floe_group);
}


TEMPLATE_PB
void PROBLEM::solve(value_type end_time, value_type dt_default, value_type out_step, bool reset){
    if (reset) create_optim_vars();
    m_domain.set_out_step(out_step);
    m_domain.set_default_time_step(dt_default);
    output_datas(); // Initial state out
    detect_proximity(); // First proximity detection
    while (m_domain.time() < end_time)
    {   
        // auto t_start = std::chrono::high_resolution_clock::now();
        step_solve();
        // auto t_end = std::chrono::high_resolution_clock::now();
        // std::cout << "Chrono STEP : " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << std::endl;
        if (*QUIT) break; // exit normally after SIGINT
    }
    std::cout << " NB STEPS : " << m_step_nb << std::endl;
}


TEMPLATE_PB
void PROBLEM::step_solve(){
    auto t0 = std::chrono::high_resolution_clock::now();
    manage_collisions();
    auto t1 = std::chrono::high_resolution_clock::now();
    compute_time_step();
    auto t2 = std::chrono::high_resolution_clock::now();
    safe_move_floe_group();
    auto t3 = std::chrono::high_resolution_clock::now();
    std::cout << "Chrono : collisions " << std::chrono::duration<double, std::milli>(t1-t0).count() << " ms + "
    << "time_step " << std::chrono::duration<double, std::milli>(t2-t1).count() << " ms + "
    << "move " << std::chrono::duration<double, std::milli>(t3-t2).count() << " ms = "
    << std::chrono::duration<double, std::milli>(t3-t0).count() << " ms" << std::endl;

    output_datas();
    m_step_nb++;
}

TEMPLATE_PB
void PROBLEM::safe_move_floe_group(){
    m_floe_group.backup_step_states();
    move_floe_group();
    while (!m_proximity_detector.update())
    {
        std::cout << "INTER "; 
        m_floe_group.recover_previous_step_states();
        m_domain.rewind_time();
        compute_time_step(); // will only divide previous time step
        if (m_domain.time_step() < m_domain.default_time_step() / 1e8)
        {   
            // Hack to bypass repeating interpenetrations...
            m_out_manager.flush();
            recover_states_from_file(m_out_manager.out_file_name(), m_domain.time() + 1);
            std::cout << "dt too small -> RECOVER STATES FROM OUT FILE" << std::endl;
            continue;
        }
        move_floe_group();
    }    
}

TEMPLATE_PB
void PROBLEM::output_datas(){
    std::cout << "----" << std::endl;
    std::cout << " Time : " << m_domain.time();
    std::cout << " | delta_t : " << m_domain.time_step();
    std::cout << " | Kinetic energy : " << m_floe_group.kinetic_energy() << std::endl;
    
    // ouput data
    if (m_domain.need_step_output())
    {
        m_out_manager.save_step(m_domain.time(), m_dynamics_manager);
        m_domain.update_next_out_limit();
    }
}


TEMPLATE_PB
void PROBLEM::detect_proximity(){
    if (m_domain.time_step() < m_domain.default_time_step() / 1e8)
    {   
        // Hack to bypass repeating interpenetrations...
        m_out_manager.flush();
        recover_states_from_file(m_out_manager.out_file_name(), m_domain.time() + 1);
        std::cout << "dt too small -> RECOVER STATES FROM OUT FILE" << std::endl;
    }
    if (!m_proximity_detector.update()) // we have a floe interpenetration
        m_domain.rewind_time();
}


TEMPLATE_PB
int PROBLEM::manage_collisions(){
    int nb_lcp = m_collision_manager.solve_contacts(m_proximity_detector.contact_graph());
    m_proximity_detector.clean_dist_opt();
    return nb_lcp;
}


TEMPLATE_PB
void PROBLEM::compute_time_step(){
    m_time_scale_manager.delta_t_secu(&m_domain);
}


TEMPLATE_PB
void PROBLEM::move_floe_group(){
    m_dynamics_manager.move_floes(m_floe_group, m_domain.time_step());
    m_domain.update_time();
}

TEMPLATE_PB
void PROBLEM::make_input_file(){
    m_out_manager.make_input_file(m_dynamics_manager);
}


}} // namespace floe::problem


#endif // PROBLEM_PROBLEM_HPP