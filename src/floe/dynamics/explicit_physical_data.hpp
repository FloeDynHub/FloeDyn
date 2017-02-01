/*!
 * \file dynamics/explicit_physical_data.hpp
 * \brief Physical datas (air & ocean speed)
 * \author Quentin Jouet
 */

#ifndef OPE_EXPLICIT_PHYSICAL_DATA_HPP
#define OPE_EXPLICIT_PHYSICAL_DATA_HPP

#include "floe/geometry/arithmetic/arithmetic.hpp"
#include "floe/geometry/arithmetic/point_operators.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
#include <ctime>
// #include <iostream> // DEBUG


namespace floe { namespace dynamics
{

/*! PhysicalData
 *
 * Gives air and ocean datas throughout time and space
 *
 */


template <typename TPoint>
class ExplicitPhysicalData
{

public:

    using point_type = TPoint;
    using value_type = decltype(TPoint::x);
    using point_vector = std::vector<point_type>;

    //! Constructor
    ExplicitPhysicalData(value_type const& time_ref) : m_time_ref{time_ref},
        m_window_width{1}, m_window_height{1},
        m_vortex_radius{0}, m_vortex_origin{0,0}, m_vortex_speed{0,0}, m_vortex_max_norm{0},
        m_water_mode{0}, m_air_mode{4} {}

    point_type water_speed(point_type pt = {0,0}) {
        return get_speed(pt, m_water_mode);
    }

    point_type air_speed(point_type pt = {0,0}) {
        return get_speed(pt, m_air_mode);
    }

    void update_water_speed(point_type diff_speed){}
    //! OBL speed accessor for output
    inline point_type OBL_speed() const { return {0,0}; }

    //! Load ocean and wind data from a topaz file
    void load_matlab_topaz_data(std::string const& filename) {}

    void set_window_size(value_type width, value_type height) {
        m_window_width = width;
        m_window_height = height;
    }

    void set_modes(int water_mode, int air_mode) {
        m_water_mode = water_mode;
        m_air_mode = air_mode;
    }

    void set_storm_mode(){
        init_random_vortex();
        m_air_mode = 5;
        m_water_mode = 0;
    }

    void init_random_vortex(){
        value_type avg_Rc{300000}, delta_Rc{110000};
        value_type avg_V{8}, delta_V{4};
        value_type min_Umax{15}, max_Umax{30};
        auto dist_Rc = std::uniform_real_distribution<double>{
            avg_Rc - delta_Rc, avg_Rc + delta_Rc};
        auto dist_V = std::uniform_real_distribution<double>{
            avg_V - delta_V, avg_V + delta_V};
        auto dist_Umax = std::uniform_real_distribution<double>{
            min_Umax, max_Umax};
        auto dist_angle = std::uniform_real_distribution<double>{
            0, 2 * M_PI};
        auto gen = std::default_random_engine{};
        gen.seed(std::time(0));
        m_vortex_radius = dist_Rc(gen);
        m_vortex_max_norm = dist_Umax(gen);
        value_type theta = dist_angle(gen);
        // 500km from {0,0} todo : less rigid choice !
        m_vortex_origin = 500000 * point_type{cos(theta), sin(theta)};
        m_vortex_speed = - dist_V(gen) * point_type{cos(theta), sin(theta)};
    }

private:

    value_type const& m_time_ref; //!< reference to time variable
    // window dimension (for generator)
    value_type m_window_width;
    value_type m_window_height;
    // vortex attributes
    value_type m_vortex_radius;
    point_type m_vortex_origin;
    point_type m_vortex_speed;
    value_type m_vortex_max_norm;
    // modes
    int m_water_mode;
    int m_air_mode;

    point_type get_speed(point_type pt = {0,0}, int mode=0){
        switch(mode){
            case 1:
                return centered_convergent_field(pt);
            case 2:
                return convergent_outside_window_field(pt);
            case 3:
                return convergent_outside_window_circular_inside_field(pt);
            case 4:
                return x_convergent_then_constant(pt);
            case 5:
                return vortex(pt);
            default:
                return {0,0};
        }
    }

    //! center convergent current (negative coeff will give divergent field)
    point_type centered_convergent_field(point_type pt = {0,0}, value_type coeff = 1) { return - coeff * pt / norm2(pt); }

    //! convergent current outside null rectangle window
    point_type convergent_outside_window_field(point_type pt = {0,0}, value_type speed = 30) {
        value_type x{0}, y{0};
        if (std::abs(pt.x) > m_window_width / 2) x = - speed * pt.x / std::abs(pt.x);
        if (std::abs(pt.y) > m_window_height / 2) y = - speed * pt.y / std::abs(pt.y);
        return {x,y};
    }

    //! convergent current outside null rectangle window
    point_type circular_inside_window_field(point_type pt = {0,0}, value_type in_speed=10) {
        value_type x{0}, y{0};
        if (std::abs(pt.x) <= m_window_width / 2 and std::abs(pt.y) <= m_window_height / 2){
            if (std::abs(pt.x) < std::abs(pt.y)) x = - 10 * pt.y / std::abs(pt.y);
            if (std::abs(pt.x) >= std::abs(pt.y)) y = 10 * pt.x / std::abs(pt.x);
        }
        return {x,y};
    }

    //! convergent current outside null rectangle window
    point_type convergent_outside_window_circular_inside_field(
        point_type pt = {0,0}, value_type out_speed=30, value_type in_speed=10
    ){
        return convergent_outside_window_field(pt, out_speed)
            + circular_inside_window_field(pt, in_speed);
    }

    //! x axis convergent wind, then constant
    point_type x_convergent_then_constant(point_type pt = {0,0}, value_type coeff = 10, point_type constant = {10,0}) {
        if (m_time_ref < 1500)
            return {0, - coeff * (pt.y / (100 + std::abs(pt.y)))};
        else
            return constant;
    }

    point_type vortex_center(){
        return m_vortex_origin + m_time_ref * m_vortex_speed;
    }

    //! x axis convergent wind, then constant
    point_type vortex(point_type pt) {
        point_type vortex_center_to_pt = pt - vortex_center();
        value_type distance_to_center = norm2(vortex_center_to_pt);
        if (distance_to_center < m_vortex_radius){
            return (m_vortex_max_norm / m_vortex_radius) *
                floe::geometry::direct_orthogonal(vortex_center_to_pt);
        } else {
            return m_vortex_max_norm * (m_vortex_radius / distance_to_center) *
                floe::geometry::direct_orthogonal(vortex_center_to_pt) / distance_to_center;
        }
    }
};


}} // namespace floe::dynamics


#endif // OPE_EXPLICIT_PHYSICAL_DATA_HPP