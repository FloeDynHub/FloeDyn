/*!
 * \file dynamics/explicit_physical_data.hpp
 * \brief Physical datas (air & ocean speed)
 * \author Quentin Jouet
 */

#ifndef OPE_EXPLICIT_PHYSICAL_DATA_HPP
#define OPE_EXPLICIT_PHYSICAL_DATA_HPP

#include "floe/geometry/arithmetic/point_operators.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
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
        m_window_width{1}, m_window_height{1}, m_water_mode{2}, m_air_mode{0} {}

    point_type water_speed(point_type pt = {0,0}) {
        switch(m_water_mode){
            case 1:
                return centered_convergent_field(pt, 0.1);
            case 2:
                return convergent_outside_window_field(pt);
            default:
                return {0,0};
        }
    }

    point_type air_speed(point_type pt = {0,0}) {
        switch(m_air_mode){
            case 1:
                return x_convergent_then_constant(pt);
            case 2:
                return convergent_outside_window_field(pt);
            default:
                return {0,0};
        }
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

private:

    value_type const& m_time_ref; //!< reference to time variable
    value_type m_window_width;
    value_type m_window_height;
    int m_water_mode;
    int m_air_mode;

    //! center convergent current (negative coeff will give divergent field)
    point_type centered_convergent_field(point_type pt = {0,0}, value_type coeff = 1) { return - coeff * pt / norm2(pt); }

    //! convergent current outside null rectangle window
    point_type convergent_outside_window_field(point_type pt = {0,0}) {
        value_type x{0}, y{0};
        if (std::abs(pt.x) > m_window_width / 2) x = - pt.x / std::abs(pt.x);
        if (std::abs(pt.y) > m_window_height / 2) y = - pt.y / std::abs(pt.y);
        return {x,y};
    }

    //! x axis convergent wind, then constant
    point_type x_convergent_then_constant(point_type pt = {0,0}, value_type coeff = 10, point_type constant = {10,0}) {
        if (m_time_ref < 1500)
            return {0, - coeff * (pt.y / (100 + std::abs(pt.y)))};
        else
            return constant;
    }
};


}} // namespace floe::dynamics


#endif // OPE_EXPLICIT_PHYSICAL_DATA_HPP