/*!
 * \file ope/external_forces.hpp
 * \brief external forces
 * \author Quentin Jouet
 */


#ifndef OPE_EXTERNAL_FORCES_HPP
#define OPE_EXTERNAL_FORCES_HPP

#include "floe/geometry/arithmetic/arithmetic.hpp"
#include "floe/geometry/arithmetic/point_operators.hpp"
#include "floe/ope/physical_data.hpp"
#include <cmath>


namespace floe { namespace ope
{

/*! ExternalForces
 *
 * Provides lambda function representing external forces for a floe
 *
 */

namespace fg = floe::geometry;

template <typename TFloeGroup>
class ExternalForces
{

public:
    using floe_group_type = TFloeGroup;
    using floe_type = typename floe_group_type::floe_type;
    using point_type = typename floe_type::point_type;
    using value_type = typename floe_type::value_type;
    using physical_data_type = PhysicalData<point_type>;

    ExternalForces(value_type const& time_ref) : m_physical_data{time_ref} {}

    std::function<point_type (value_type, value_type)> total_drag(floe_type& floe);
    std::function<value_type (value_type, value_type)> total_rot_drag(floe_type& floe);
    point_type coriolis_effect(floe_type& floe);

    inline void load_matlab_topaz_data(std::string const& filename) {
        m_physical_data.load_matlab_topaz_data(filename);
    }

    inline value_type OBL_surface_mass() const { return h_w * rho_w; }
    inline void update_water_speed( point_type diff_speed ) { m_physical_data.update_water_speed( diff_speed ); }

    // forces applied on floes
    std::function<point_type (value_type, value_type)> ocean_drag_2(floe_type& floe);

    // forces applied on ocean
    point_type air_drag_ocean();
    point_type ocean_coriolis(point_type p);
    point_type deep_ocean_friction();

private:

    physical_data_type m_physical_data;

    const value_type rho_w = 1024.071; // Water density
    const value_type C_w = 5 * 1e-3; // Oceanic skin drag coeff
    const value_type rho_a = 1.341; // Air density
    const value_type C_a = 1.7 * 1e-3; // Atmospheric skin drag coeff

    const value_type R_earth = 6371 * 1e3; // earth radius
    const value_type V_earth = 7.292 * 1e-5; // Earth angular velocity

    const value_type O_latitude = 80.207; // Origin latitude

    const value_type gamma = 1e-5; // m/s friction velocity within the OBL
    const value_type h_w = 15; // OBL height

    value_type coriolis_coeff(point_type p);

    inline point_type water_speed(point_type p){ return m_physical_data.water_speed(p); }
    inline point_type air_speed(point_type p){ return m_physical_data.air_speed(p); }

    // forces applied on floes
    std::function<point_type (point_type&)> ocean_drag(floe_type& floe);
    std::function<point_type (point_type&)> air_drag();

};


template<typename TFloeGroup>
using value = typename TFloeGroup::floe_type::value_type;


template <typename TFloeGroup>
std::function<typename TFloeGroup::floe_type::point_type (
    typename TFloeGroup::floe_type::point_type&)>
ExternalForces<TFloeGroup>::ocean_drag(floe_type& floe)
{   
    auto& state = floe.state();
    return [&](point_type& p)
    {
        auto speed_p = state.speed 
            + state.rot * fg::direct_orthogonal(p - state.pos);
        auto V = water_speed(p) - speed_p;
        return rho_w * C_w * norm2(V) * V;
    };
}

template <typename TFloeGroup>
std::function<typename TFloeGroup::floe_type::point_type (
    value<TFloeGroup>, value<TFloeGroup>)>
ExternalForces<TFloeGroup>::ocean_drag_2(floe_type& floe)
{   
    return [&](value_type x, value_type y)
    {
        point_type p{x,y};
        return ocean_drag(floe)(p);
    };
}


template <typename TFloeGroup>
std::function<typename TFloeGroup::floe_type::point_type (
    typename TFloeGroup::floe_type::point_type&)>
ExternalForces<TFloeGroup>::air_drag()
{
    return [&](point_type& p)
    {
        auto f = air_speed(p);
        return rho_a * C_a * norm2(f) * f;
    };
}


template <typename TFloeGroup>
typename TFloeGroup::floe_type::point_type
ExternalForces<TFloeGroup>::coriolis_effect(floe_type& floe)
{
    return - coriolis_coeff(floe.state().pos) * fg::direct_orthogonal(floe.state().speed);
}

template <typename TFloeGroup>
value<TFloeGroup>
ExternalForces<TFloeGroup>::coriolis_coeff(point_type p)
{
    auto phi = (O_latitude * M_PI / 180 + p.y / R_earth); // in radian
    return 2 * V_earth * sin(phi);
}


template <typename TFloeGroup>
std::function<typename TFloeGroup::floe_type::point_type (
    value<TFloeGroup>, value<TFloeGroup>)>
ExternalForces<TFloeGroup>::total_drag(floe_type& floe)
{
    return [&](value_type x, value_type y)
    {
        point_type p{x,y};
        return ocean_drag(floe)(p) + air_drag()(p);
    };
}


template <typename TFloeGroup>
std::function<value<TFloeGroup> (
    value<TFloeGroup>, value<TFloeGroup>)>
ExternalForces<TFloeGroup>::total_rot_drag(floe_type& floe)
{
    return [&](value_type x, value_type y)
    {
        point_type p{x,y};
        return fg::cross_product_value(p - floe.state().pos, total_drag(floe)(x, y));
    };
}


template <typename TFloeGroup>
typename ExternalForces<TFloeGroup>::point_type
ExternalForces<TFloeGroup>::air_drag_ocean()
{   
    auto f = air_speed({0,0}); // Wind is uniform in space for now
    return rho_a * C_a * norm2(f) * f;
}

template <typename TFloeGroup>
typename ExternalForces<TFloeGroup>::point_type
ExternalForces<TFloeGroup>::ocean_coriolis(point_type p)
{   
    return - coriolis_coeff(p) * fg::direct_orthogonal(m_physical_data.water_speed());
}

template <typename TFloeGroup>
typename ExternalForces<TFloeGroup>::point_type
ExternalForces<TFloeGroup>::deep_ocean_friction()
{   
    return - ( gamma / h_w ) * m_physical_data.water_speed();
}


}} // namespace floe::ope


#endif // OPE_EXTERNAL_FORCES_HPP