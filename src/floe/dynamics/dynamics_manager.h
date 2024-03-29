/*!
 * \file dynamics/dynamics_manager.h
 * \brief Dynamics manager
 * \author Quentin Jouet
 */

#ifndef OPE_DYNAMICS_MANAGER_H
#define OPE_DYNAMICS_MANAGER_H

#ifdef _OPENMP
#include <omp.h>
#endif

 #include <iostream> // DEBUG
#include <random>


namespace floe { namespace dynamics
{

/*! DynamicsManager
 *
 * Operator for dynamics processing (Floes and ocean kinematics)
 *
 */


template <typename TExternalForces, typename TFloeGroup>
class DynamicsManager
{

public:
    using external_forces_type = TExternalForces;
    using floe_group_type = TFloeGroup;
    using floe_type = typename floe_group_type::floe_type;
    using point_type = typename floe_type::point_type;
    using real_type = typename floe_type::real_type;
    using state_type = typename floe_type::state_type;

    //! Constructor
    DynamicsManager(real_type const& time_ref, int OBL_status) : m_external_forces{time_ref}, m_ocean_window_area{0},
        m_OBL_status{OBL_status}, m_rand_speed_add{false}, m_rand_norm{1e-7} {}

    //! Floes state update (+ update_ocean call, returns update_ocean's return value)
    point_type move_floes(floe_group_type& floe_group, real_type delta_t);
    //! Ocean state update, returns difference speed applied
    point_type update_ocean(floe_group_type& floe_group, real_type delta_t, point_type floes_force = {0,0});

    //! Load ocean and wind data from a topaz file
    inline void load_matlab_topaz_data(std::string const& filename) {
        m_external_forces.load_matlab_topaz_data(filename);
    }

    //! Load ocean window area (box surrounding floes)
    void load_matlab_ocean_window_data(std::string const& filename, floe_group_type const& floe_group);

    //! for output
    inline point_type OBL_speed() const { return m_external_forces.OBL_speed(); }
    inline void set_OBL_speed(point_type OBL_speed) { 
        return this->m_external_forces.get_physical_data().set_OBL_speed(OBL_speed); }
    inline void set_OBL_status(int status) { m_OBL_status = status; }
    //! Ocean window area setter
    inline void set_ocean_window_area(real_type area) { m_ocean_window_area = area; }

    //!< extra random velocities
    inline void set_rand_speed_add(bool rand_speed_add) {m_rand_speed_add = rand_speed_add;}
    inline void set_norm_rand_speed(real_type rand_norm) {m_rand_norm = rand_norm;}

    //! Accessor for specific use
    external_forces_type& get_external_forces() { return m_external_forces; }

protected:

    external_forces_type m_external_forces; //! External forces manager
    real_type m_ocean_window_area; //! Ocean window area (for OBL computing)
    int m_OBL_status; //! OBL (Oceanic Boundary Layer) status : 0 = no coupling, 1 = coupling
    std::default_random_engine m_random_generator;

    bool m_rand_speed_add; //!< extra random velocities 
    real_type m_rand_norm; //!< norm of these extra random velocities

    //! Move one floe
    virtual void move_floe(floe_type& floe, real_type delta_t);
    //! Ocean window area accessor
    virtual real_type ocean_window_area() { return m_ocean_window_area; }
};


}} // namespace floe::dynamics


#endif // OPE_DYNAMICS_MANAGER_H