/*!
 * \file problem/problem.hpp
 * \brief Smooth problem
 * \author Quentin Jouet
 */


#ifndef PROBLEM_PROBLEM_HPP
#define PROBLEM_PROBLEM_HPP

#include "floe/io/hdf5_manager.h"
// #include "floe/io/false_hdf5_manager.hpp"
#include "floe/io/multi_out_manager.hpp"

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
    #ifdef MULTIOUTPUT
        using out_manager_type = io::MultiOutManager<io::HDF5Manager<TFloeGroup, TDynamicsManager>>;
    #else
        using out_manager_type = io::HDF5Manager<TFloeGroup, TDynamicsManager>;
    #endif

    using real_type = typename TFloeGroup::real_type;
    using point_type = typename TFloeGroup::point_type;
    using floe_group_type = TFloeGroup; // generator accessor
    using time_scale_manager_type = domain::TimeScaleManager<typename TProxymityDetector::proximity_data_type>;
    using proximity_detector_type = TProxymityDetector;

    //! Default constructor.
    Problem(real_type epsilon=0.4, int OBL_status=0);

    //! Solver of the problem (main method)
    virtual void solve(real_type end_time, real_type dt_default, real_type out_step = 0, bool reset = true);

    virtual void load_config(std::string const& filename);
    //! Load ocean and wind data from a topaz file
    inline void load_matlab_topaz_data(std::string const& filename) {
        m_dynamics_manager.load_matlab_topaz_data(filename);
    }
    //! Recover simulation state from previous ouput file, at any recorded time t
    virtual void recover_states_from_file(std::string const& filename, real_type t, bool keep_as_outfile=true);

    //! set existing floe_group (from generator for example)
    virtual void set_floe_group(floe_group_type& floe_group);
    //! Floe group accessor for config generator
    inline TFloeGroup& get_floe_group() { return m_floe_group; }
    //! Dynamics mgr accessor for config generator
    inline TDynamicsManager& get_dynamics_manager(){ return m_dynamics_manager; }
    
    /*! Floe Concentration
     *
     * \return the floe concentration using FloeGroup::floe_concentration. (\f$ 0 <= \f$ floe concentration \f$ <= 1 \f$).
     */
    virtual real_type floe_concentration() { return m_floe_group.floe_concentration(); }
    void make_input_file();
    //! Access detector
    inline proximity_detector_type& proximity_detector() { return m_proximity_detector; }
    //! Initializing proximity detector with floe set
    void create_optim_vars();
    //!< OutManager accessor
    inline out_manager_type& get_out_manager() {return m_out_manager;}
    //!< LCP solver accessor
    inline TCollisionManager& get_lcp_manager() { return m_collision_manager; }

    const std::atomic<bool>* QUIT; //!< Exit signal

    // to be used in generation mode, to allow different behaviour
    inline void set_is_generator() {m_is_generator = true;};


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
    point_type move_floe_group();
    //! Handle output_datas (console + out file)
    void output_datas();

    // to allow different behaviour in the generation phase
    bool m_is_generator;
    // run time breakdown  
    std::chrono::duration<double, std::nano> m_collisionTime;
    std::chrono::duration<double, std::nano> m_timeStepTime;
    std::chrono::duration<double, std::nano> m_moveTime;
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
PROBLEM::Problem(real_type epsilon, int OBL_status) :
        QUIT{nullptr},
        m_domain{},
        m_proximity_detector{},
        m_collision_manager{epsilon},
        m_dynamics_manager{m_domain.time(), OBL_status},
        m_floe_group{},
        m_step_nb{0},
        m_out_manager{m_floe_group},
        m_is_generator{false},
        m_collisionTime{},
        m_timeStepTime{},
        m_moveTime{}
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
void PROBLEM::recover_states_from_file(std::string const& filename, real_type t, bool keep_as_outfile){
    real_type saved_time = m_out_manager.recover_states(filename, t, m_floe_group, m_dynamics_manager, keep_as_outfile);
    std::cout << "RECOVER : " << saved_time << std::endl;
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
    m_proximity_detector.set_floe_group(m_floe_group);
}


TEMPLATE_PB
void PROBLEM::solve(real_type end_time, real_type dt_default, real_type out_step, bool reset){
    if (reset) this->create_optim_vars();
    this->m_domain.set_default_time_step(dt_default);
    this->m_out_manager.set_out_step(out_step, this->m_domain.time());
    this->output_datas(); // Initial state out
    this->detect_proximity(); // First proximity detection
    while (this->m_domain.time() < end_time)
    {   
        // auto t_start = std::chrono::high_resolution_clock::now();
        this->step_solve();
        // auto t_end = std::chrono::high_resolution_clock::now();
        // std::cout << "Chrono STEP : " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << std::endl;
        if (*this->QUIT) break; // exit normally after SIGINT
    }
    std::cout << " NB STEPS : " << this->m_step_nb << std::endl;
    double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(m_collisionTime).count();
    std::cout << " total collision time : " << time_taken*1e-9 << " s" << std::endl;
    std::cout << " total time step time : " << std::chrono::duration_cast<std::chrono::nanoseconds>(m_timeStepTime).count()*1e-9 << " s" << std::endl;
    std::cout << " total move time : " << std::chrono::duration_cast<std::chrono::nanoseconds>(m_moveTime).count()*1e-9 << " s" << std::endl;
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
    if (this->m_dynamics_manager.get_external_forces().get_physical_data().get_air_mode()==5 || this->m_dynamics_manager.get_external_forces().get_physical_data().get_air_mode()==6) { //!< only if the external forces is a vortex
        for (size_t i=0; i<m_dynamics_manager.get_external_forces().get_physical_data().get_nb_vortex(); ++i) {
            std::cout << "the vortex wind speed is: " << 
                this->m_dynamics_manager.get_external_forces().get_physical_data().get_vortex_wind_speed(i) 
                << std::endl;
        }
        /*
        std::vector<real_type> statsKinematic = m_floe_group.get_KinematicFloeWithMaxKineticEnergy();
        std::cout << "Velocities from Floe Id: " << statsKinematic[0] <<
        " At position x= " << statsKinematic[1] << ", y= " << statsKinematic[2] <<
        " with the max Kinetic Energy: Vx = " <<
        statsKinematic[3] << ", Vy = " << statsKinematic[4] << ", Vtheta = " <<
        statsKinematic[5] << std::endl;

        std::vector<real_type> statsDist = m_proximity_detector.get_statsDistance();
        std::cout << "\nMinimal true distance = " << statsDist[0] << " between floes: ("
        << statsDist[1] << ", " << statsDist[2] << ")" << "far as: " << statsDist[3]
        << std::endl;
        std::cout << "optimal distance = " << statsDist[4] << " between floes: ("
        << statsDist[5] << ", " << statsDist[6] << ")" << "far as: " << statsDist[7]
        << std::endl;
        std::vector<real_type> statsDeltaT = m_time_scale_manager.get_statsDeltaT();
        std::cout << "\nMinimal Fast Delta T = " << statsDeltaT[0] << " between floes: ("
        << statsDeltaT[1] << ", " << statsDeltaT[2] << ")" << "far as: " << statsDeltaT[3]
        << std::endl;
        std::cout << "Minimal Safe Delta T = " << statsDeltaT[4] << " between floes: ("
        << statsDeltaT[5] << ", " << statsDeltaT[6] << ")" << "far as: " << statsDeltaT[7]
        << std::endl;
         */
    }
    std::cout << "Chrono : collisions " << std::chrono::duration<double, std::milli>(t1-t0).count() << " ms + "
    << "time_step " << std::chrono::duration<double, std::milli>(t2-t1).count() << " ms + "
    << "move " << std::chrono::duration<double, std::milli>(t3-t2).count() << " ms = "
    << std::chrono::duration<double, std::milli>(t3-t0).count() << " ms" << std::endl;

    m_collisionTime+=std::chrono::duration<double, std::nano>(t1-t0);
    m_timeStepTime+=std::chrono::duration<double, std::nano>(t2-t1);
    m_moveTime+=std::chrono::duration<double, std::nano>(t3-t2);

    output_datas();
    m_step_nb++;
}

TEMPLATE_PB
void PROBLEM::safe_move_floe_group(){
    m_floe_group.backup_step_states();
    move_floe_group();
    real_type norm_rand_speed(0);
    while (!m_proximity_detector.update())
    {
        std::cout << "INTER "; 
        m_floe_group.recover_previous_step_states();
        m_domain.rewind_time();
        compute_time_step(); // will only divide previous time step
        m_dynamics_manager.set_norm_rand_speed(1e-7); 
        if (m_domain.time_step() < m_domain.default_time_step() / 1e5) // 1e8 from Q.Jouet
        {   
            // Hack to bypass repeating interpenetrations...
            m_out_manager.flush();
            recover_states_from_file(m_out_manager.out_file_name(), m_domain.time() + 1);
            std::cout << "dt too small -> RECOVER STATES FROM OUT FILE" << std::endl;

            // adding a random velocity component to avoid infinite loops
            // in generator mode the value is between 0 and 1. 
            // in physical run mode the value is chosen smaller in order not to impair the physical sense of the solution :
            // this epsilon is chosen depending on the mean velocity, the maximum floe size and the current time step.   
            if (m_is_generator) 
            {
                norm_rand_speed = static_cast <real_type> (std::rand()) / static_cast <real_type> (RAND_MAX);
                std::cout << "Adding a random component to velocities of " << norm_rand_speed << std::endl;
                m_dynamics_manager.set_rand_speed_add(true);
                m_dynamics_manager.set_norm_rand_speed(norm_rand_speed); 
            }
            // else
            // {
            //     std::cout << "Adding a random epsilon to velocities to all floes " << std::endl;
            //     norm_rand_speed = static_cast <real_type> (std::rand()) / static_cast <real_type> (RAND_MAX);
            //     m_dynamics_manager.set_rand_speed_add(true);
            //     m_dynamics_manager.set_norm_rand_speed(norm_rand_speed); 
            // } 
            continue;
        }
        move_floe_group();
        m_dynamics_manager.set_norm_rand_speed(1e-7); 
        // m_dynamics_manager.set_rand_speed_add(false); // it seems SimuRunner sets it to true by default 
    }    
}

TEMPLATE_PB
void PROBLEM::output_datas(){
    std::cout << "----" << std::endl;
    std::cout << " Time : " << this->m_domain.time();
    std::cout << " | delta_t : " << this->m_domain.time_step();
    std::cout << " | Kinetic energy : " << this->m_floe_group.kinetic_energy() << std::endl;
    // ouput data
    m_out_manager.save_step_if_needed(this->m_domain.time(), this->m_dynamics_manager);
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
typename TFloeGroup::point_type PROBLEM::move_floe_group(){
    point_type resp = m_dynamics_manager.move_floes(m_floe_group, m_domain.time_step());
    m_domain.update_time();
    return resp;
}

TEMPLATE_PB
void PROBLEM::make_input_file(){
    m_out_manager.make_input_file(m_dynamics_manager);
}

}} // namespace floe::problem


#endif // PROBLEM_PROBLEM_HPP