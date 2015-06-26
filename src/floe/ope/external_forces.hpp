/*!
 * \file ope/external_forces.hpp
 * \brief external forces
 * \author Quentin Jouet
 */

 //TODO : better calcul repartition
 //       create class for physical datas (air and wind)

#ifndef OPE_EXTERNAL_FORCES_HPP
#define OPE_EXTERNAL_FORCES_HPP

#include "floe/geometry/arithmetic/arithmetic.hpp"
#include "floe/geometry/arithmetic/point_operators.hpp"
#include <math.h>


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

    std::function<point_type (value_type, value_type)> total_drag(floe_type& floe);
    std::function<point_type (value_type, value_type)> total_drag2(floe_type& floe);
    std::function<value_type (value_type, value_type)> total_rot_drag(floe_type& floe);
    point_type coriolis_effect(floe_type& floe);

private:


    // PROVISOIRE
    inline point_type water_speed(point_type){ return point_type{0,1}; }
    inline point_type air_speed(point_type){ return point_type{-0,0}; }
    // PROVISOIRE

    const value_type rho_w = 1024.071; // Water density
    const value_type C_w = 5 * 1e-3; // Oceanic skin drag coeff
    const value_type rho_a = 1.341; // Air density
    const value_type C_a = 1.7 * 1e-3; // Atmospheric skin drag coeff

    const value_type R_earth = 6371 * 1e3; // earth radius
    const value_type V_earth = 7.292 * 1e-5; // Earth angular velocity

    const value_type O_latitude = 80.207; // Origin latitude

    value_type coriolis_coeff(floe_type& floe);

    void move_floe(floe_type& floe, value_type delta_t);
    std::function<point_type (point_type&)> ocean_drag(floe_type& floe);
    point_type ocean_drag2(floe_type& floe, point_type& p);
    std::function<point_type (point_type&)> air_drag(floe_type& floe);
    point_type air_drag2(floe_type& floe, point_type& p);

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
typename TFloeGroup::floe_type::point_type
ExternalForces<TFloeGroup>::ocean_drag2(floe_type& floe, point_type& p)
{
    auto& state = floe.state();
    auto V = water_speed(p) - state.speed - state.rot * fg::direct_orthogonal(p - state.pos);
    return rho_w * C_w * norm2(V) * V;
}


template <typename TFloeGroup>
std::function<typename TFloeGroup::floe_type::point_type (
    typename TFloeGroup::floe_type::point_type&)>
ExternalForces<TFloeGroup>::air_drag(floe_type& floe)
{
    return [&](point_type& p)
    {
        auto f = air_speed(p);
        return rho_a * C_a * norm2(f) * f;
    };
}


template <typename TFloeGroup>
typename TFloeGroup::floe_type::point_type
ExternalForces<TFloeGroup>::air_drag2(floe_type& floe, point_type& p)
{
    auto f = air_speed(p);
    return rho_a * C_a * norm2(f) * f;
}


template <typename TFloeGroup>
typename TFloeGroup::floe_type::point_type
ExternalForces<TFloeGroup>::coriolis_effect(floe_type& floe)
{
    return - coriolis_coeff(floe) * fg::direct_orthogonal(floe.state().speed);
}

template <typename TFloeGroup>
value<TFloeGroup>
ExternalForces<TFloeGroup>::coriolis_coeff(floe_type& floe)
{
    auto phi = (O_latitude * M_PI / 180 + floe.state().pos.y / R_earth); // in radian
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
        return ocean_drag(floe)(p) + air_drag(floe)(p);
    };
}

template <typename TFloeGroup>
std::function<typename TFloeGroup::floe_type::point_type (
    value<TFloeGroup>, value<TFloeGroup>)>
ExternalForces<TFloeGroup>::total_drag2(floe_type& floe)
{
    return [&](value_type x, value_type y)
    {
        point_type p{x,y};
        return ocean_drag2(floe, p) + air_drag2(floe, p);
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
        return fg::cross_product_value(p - floe.state().pos, total_drag2(floe)(x, y));
    };
}


}} // namespace floe::ope


#endif // OPE_EXTERNAL_FORCES_HPP