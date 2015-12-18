/*!
 * \file ope/dynamics_manager.h
 * \brief Dynamics manager
 * \author Quentin Jouet
 */

#ifndef OPE_DYNAMICS_MANAGER_H
#define OPE_DYNAMICS_MANAGER_H

#ifdef _OPENMP
#include <omp.h>
#endif

 #include <iostream> // DEBUG


namespace floe { namespace ope
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
    using value_type = typename floe_type::value_type;
    using state_type = typename floe_type::state_type;

    //! Constructor
    DynamicsManager(value_type const& time_ref) : m_external_forces{time_ref}, m_ocean_window_area{0}, m_OBL_status{0} {}

    //! Floes state update
    void move_floes(floe_group_type& floe_group, value_type delta_t);
    //! Ocean state update
    void update_ocean(floe_group_type& floe_group, value_type delta_t);

    //! Load ocean and wind data from a topaz file
    inline void load_matlab_topaz_data(std::string const& filename) {
        m_external_forces.load_matlab_topaz_data(filename);
    }

    //! Load ocean window area (box surrounding floes)
    void load_matlab_ocean_window_data(std::string const& filename, floe_group_type const& floe_group);

    //! for output
    inline point_type OBL_speed() const { return m_external_forces.OBL_speed(); }
    inline void set_OBL_speed(point_type OBL_speed) { return m_external_forces.update_water_speed(OBL_speed); }
    inline void set_OBL_status(int status) { m_OBL_status = status; }
    //! Ocean window area setter
    inline void set_ocean_window_area(value_type area) { m_ocean_window_area = area; }

    //! Accessor for specific use
    external_forces_type& get_external_forces() { return m_external_forces; }

protected:

    external_forces_type m_external_forces; //! External forces manager
    value_type m_ocean_window_area; //! Ocean window area (for OBL computing)
    int m_OBL_status; //! OBL (Oceanic Boundary Layer) status : 0 = no coupling, 1 = coupling

    //! Move one floe
    virtual void move_floe(floe_type& floe, value_type delta_t);
    //! Ocean window area accessor
    virtual value_type ocean_window_area() { return m_ocean_window_area; }
};


}} // namespace floe::ope


#endif // OPE_DYNAMICS_MANAGER_H