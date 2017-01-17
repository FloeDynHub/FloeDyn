/*!
 * \file dynamics/external_forces.hpp
 * \brief External forces
 * \author Quentin Jouet
 */


#ifndef OPE_EXTERNAL_FORCES_HPP
#define OPE_EXTERNAL_FORCES_HPP

#include "floe/geometry/arithmetic/arithmetic.hpp"
#include "floe/geometry/arithmetic/point_operators.hpp"
#include <cmath>


namespace floe { namespace dynamics
{

/*! ExternalForces
 *
 * Provides lambda function representing external forces for a floe
 *
 */

namespace fg = floe::geometry;

template <typename TFloe, typename TPhysicalData>
class ExternalForces
{

public:
    using floe_type = TFloe;
    using point_type = typename floe_type::point_type;
    using value_type = typename floe_type::value_type;
    using physical_data_type = TPhysicalData;
    // using physical_data_type = PhysicalData<point_type>;

    ExternalForces(value_type const& time_ref) : m_physical_data{time_ref} {}

    //! Sum of different drag effects on a floe
    std::function<point_type (value_type, value_type)> total_drag(floe_type& floe);
    //! Sum of different rotational drag effects on a floe
    std::function<value_type (value_type, value_type)> total_rot_drag(floe_type& floe);
    //! Coriolis effect on a floe
    point_type coriolis_effect(floe_type& floe);

    //! Load ocean and wind data from a topaz file
    inline void load_matlab_topaz_data(std::string const& filename) {
        m_physical_data.load_matlab_topaz_data(filename);
    }

    //! Surface mass of Oceaninc Boundary Layer
    inline value_type OBL_surface_mass() const { return h_w * rho_w; }
    //! Update water speed
    inline void update_water_speed( point_type diff_speed ) { m_physical_data.update_water_speed( diff_speed ); }
    //! OBL speed accessor for output
    inline point_type OBL_speed() const { return m_physical_data.OBL_speed(); }

    // FLOES
    //! Signature Variation of ocean_drag
    std::function<point_type (value_type, value_type)> ocean_drag_2(floe_type& floe);

    // OCEAN
    //! Air drag on ocean
    point_type air_drag_ocean();
    //! Coriolis effect on ocean
    point_type ocean_coriolis(point_type p);
    //! Deep ocean friction effect on ocean
    point_type deep_ocean_friction();

    //! Accessor for specific use
    physical_data_type& get_physical_data() { return m_physical_data; }

private:

    physical_data_type m_physical_data;

    const value_type rho_w = 1024.071; //!< Water density
    const value_type C_w = 5 * 1e-3; //!< Oceanic skin drag average coeff
    const value_type rho_a = 1.341; //!< Air density
    const value_type C_a = 1.7 * 1e-3; //!< Atmospheric skin drag average coeff

    const value_type R_earth = 6371 * 1e3; //!< earth radius
    const value_type V_earth = 7.292 * 1e-5; //!< Earth angular velocity

    const value_type O_latitude = 80.207; //!< Origin latitude

    const value_type gamma = 1e-5; //!< m/s friction velocity within the OBL
    const value_type h_w = 15; //!< OBL height

    //! Coriolis coefficient at a point
    value_type coriolis_coeff(point_type p);

    //! Water speed accessor
    inline point_type water_speed(point_type p){ return m_physical_data.water_speed(p); }
    //! Air speed accessor
    inline point_type air_speed(point_type p){ return m_physical_data.air_speed(p); }

    //! Ocean drag effect on a floe
    std::function<point_type (point_type&)> ocean_drag(floe_type& floe);
    //! Air drag effect on a floe
    std::function<point_type (point_type&)> air_drag();

};


template<typename TFloe>
using value = typename TFloe::value_type;


template <typename TFloe, typename TPhysicalData>
std::function<typename TFloe::point_type (
    typename TFloe::point_type&)>
ExternalForces<TFloe, TPhysicalData>::ocean_drag(floe_type& floe)
{   
    auto& state = floe.state();
    return [&](point_type& p)
    {
        auto speed_p = state.speed 
            + state.rot * fg::direct_orthogonal(p - state.pos);
        auto V = water_speed(p) - speed_p;
        return rho_w * floe.static_floe().C_w() * norm2(V) * V;
    };
}

template <typename TFloe, typename TPhysicalData>
std::function<typename TFloe::point_type (
    value<TFloe>, value<TFloe>)>
ExternalForces<TFloe, TPhysicalData>::ocean_drag_2(floe_type& floe)
{   
    return [&](value_type x, value_type y)
    {
        point_type p{x,y};
        return ocean_drag(floe)(p);
    };
}


template <typename TFloe, typename TPhysicalData>
std::function<typename TFloe::point_type (
    typename TFloe::point_type&)>
ExternalForces<TFloe, TPhysicalData>::air_drag()
{
    return [&](point_type& p)
    {
        auto f = air_speed(p);
        return rho_a * C_a * norm2(f) * f;
    };
}


template <typename TFloe, typename TPhysicalData>
typename TFloe::point_type
ExternalForces<TFloe, TPhysicalData>::coriolis_effect(floe_type& floe)
{
    return - coriolis_coeff(floe.state().pos) * fg::direct_orthogonal(floe.state().speed);
}

template <typename TFloe, typename TPhysicalData>
value<TFloe>
ExternalForces<TFloe, TPhysicalData>::coriolis_coeff(point_type p)
{
    auto phi = (O_latitude * M_PI / 180 + p.y / R_earth); // in radian
    return 2 * V_earth * sin(phi);
}


template <typename TFloe, typename TPhysicalData>
std::function<typename TFloe::point_type (
    value<TFloe>, value<TFloe>)>
ExternalForces<TFloe, TPhysicalData>::total_drag(floe_type& floe)
{
    return [&](value_type x, value_type y)
    {
        point_type p{x,y};
        return ocean_drag(floe)(p) + air_drag()(p);
    };
}


template <typename TFloe, typename TPhysicalData>
std::function<value<TFloe> (
    value<TFloe>, value<TFloe>)>
ExternalForces<TFloe, TPhysicalData>::total_rot_drag(floe_type& floe)
{
    return [&](value_type x, value_type y)
    {
        point_type p{x,y};
        return fg::cross_product_value(p - floe.state().pos, total_drag(floe)(x, y));
    };
}


template <typename TFloe, typename TPhysicalData>
typename ExternalForces<TFloe, TPhysicalData>::point_type
ExternalForces<TFloe, TPhysicalData>::air_drag_ocean()
{   
    auto f = air_speed({0,0}); // Wind is uniform in space for now
    return rho_a * C_a * norm2(f) * f;
}

template <typename TFloe, typename TPhysicalData>
typename ExternalForces<TFloe, TPhysicalData>::point_type
ExternalForces<TFloe, TPhysicalData>::ocean_coriolis(point_type p)
{   
    return - coriolis_coeff(p) * fg::direct_orthogonal(m_physical_data.water_speed());
}

template <typename TFloe, typename TPhysicalData>
typename ExternalForces<TFloe, TPhysicalData>::point_type
ExternalForces<TFloe, TPhysicalData>::deep_ocean_friction()
{   
    return - ( gamma / h_w ) * m_physical_data.water_speed();
}


}} // namespace floe::dynamics


#endif // OPE_EXTERNAL_FORCES_HPP