/*!
 * \file dynamics/physical_data.hpp
 * \brief Physical datas (air & ocean speed)
 * \author Quentin Jouet
 */

#ifndef OPE_PHYSICAL_DATA_HPP
#define OPE_PHYSICAL_DATA_HPP

#include "floe/geometry/arithmetic/point_operators.hpp"
#include "floe/io/matlab/topaz_import.hpp"
#include <cmath>
#include <vector>
#include <iostream> // DEBUG


namespace floe { namespace dynamics
{

/*! PhysicalData
 *
 * Gives air and ocean datas throughout time and space
 *
 */


template <typename TPoint>
class PhysicalData
{

public:

    using point_type = TPoint;
    using real_type = decltype(TPoint::x);
    using point_vector = std::vector<point_type>;

    //! Constructor
    PhysicalData(real_type const& time_ref) :
        m_ocean_data_hours{}, m_air_data_hours{},
        m_ocean_data_minutes{}, m_air_data_minutes{},
        m_time_ref{time_ref},
        m_geo_relative_water_speed{0, 0},
        m_window_width{1}, m_window_height{1},
        m_vortex_radius{0}, m_vortex_origin{0,0}, m_vortex_speed{0,0}, m_vortex_max_norm{0},
        m_water_mode{0}, m_air_mode{0} {}

    //! water speed accessor
    point_type water_speed(point_type pt = {0,0});
    //! air speed accessor
    point_type air_speed(point_type pt = {0,0});
    //! OBL update speed
    void update_water_speed(point_type diff_speed);
    //! OBL speed accessor for output
    inline point_type OBL_speed() const { return m_geo_relative_water_speed; }
    void set_OBL_speed(point_type speed) { m_geo_relative_water_speed = speed; }
    //! Load ocean and wind data from a topaz file
    void load_matlab_topaz_data(std::string const& filename);
    //! For modes depending on an artificial ocean window (generator)
    void set_window_size(real_type width, real_type height) {
        m_window_width = width;
        m_window_height = height;
    }
    //! Air and water conditions mode setter
    void set_modes(int water_mode, int air_mode) {
        m_water_mode = water_mode;
        m_air_mode = air_mode;
    }
    //! Activate storm (random vortex) mode
    void set_storm_mode(){
        this->init_random_vortex();
        m_air_mode = 5;
        m_water_mode = -1;
    }
    //! Random vortex initialisation
    void init_random_vortex();

private:

    point_vector m_ocean_data_hours; //!< Geostrophic datas
    point_vector m_air_data_hours; //!< Geostrophic datas
    point_vector m_ocean_data_minutes; //!< Geostrophic datas
    point_vector m_air_data_minutes; //!< Geostrophic datas
    real_type const& m_time_ref; //!< reference to time variable

    point_type m_geo_relative_water_speed; //!< Water speed correction compared to geostrophic data

    // window dimension (for generator)
    real_type m_window_width;
    real_type m_window_height;
    // vortex attributes
    real_type m_vortex_radius;
    point_type m_vortex_origin;
    point_type m_vortex_speed;
    real_type m_vortex_max_norm;
    // modes
    int m_water_mode;
    int m_air_mode;

    //! Interpolate all hour datas to minute datas
    void interpolate_hour_to_minute();
    //! Interpolate individual hour datas to minute datas
    void interpolate_hour_to_minute(point_vector const&, point_vector&);
    //! Get value in minute datas from time
    point_type minute_value(real_type t, point_vector const&);
    //! Get geostrophic (topaz) water speed
    point_type geostrophic_water_speed(point_type = {0,0});
    //! Get topaz air speed
    point_type topaz_air_speed(point_type = {0,0});
    //! Mode router to get air or water speed
    point_type get_speed(point_type pt = {0,0}, int mode=0);
    //! center convergent current (negative coeff will give divergent field)
    point_type centered_convergent_field(point_type pt = {0,0}, real_type coeff = 1) {
        return - coeff * pt / norm2(pt);
    }
    //! convergent current outside null rectangle window
    point_type convergent_outside_window_field(point_type pt = {0,0}, real_type speed=0.1) {
        real_type x{0}, y{0};
        if (std::abs(pt.x) > m_window_width / 2) x = - speed * pt.x / std::abs(pt.x);
        if (std::abs(pt.y) > m_window_height / 2) y = - speed * pt.y / std::abs(pt.y);
        return {x,y};
    }
    //! convergent current outside null rectangle window
    point_type circular_inside_window_field(point_type pt = {0,0}, real_type in_speed=1) {
        real_type x{0}, y{0};
        if (std::abs(pt.x) <= m_window_width / 2 and std::abs(pt.y) <= m_window_height / 2){
            if (std::abs(pt.x) < std::abs(pt.y)) x = - 10 * pt.y / std::abs(pt.y);
            if (std::abs(pt.x) >= std::abs(pt.y)) y = 10 * pt.x / std::abs(pt.x);
        }
        return {x,y};
    }
    //! convergent current outside null rectangle window
    point_type convergent_outside_window_circular_inside_field(
        point_type pt = {0,0}, real_type out_speed=30, real_type in_speed=10
    ){
        return convergent_outside_window_field(pt, out_speed)
            + circular_inside_window_field(pt, in_speed);
    }
    //! x axis convergent wind, then constant
    point_type x_convergent_then_constant(point_type pt = {0,0}, real_type coeff = 10, point_type constant = {10,0}) {
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
        real_type distance_to_center = norm2(vortex_center_to_pt);
        if (distance_to_center < m_vortex_radius){
            return (m_vortex_max_norm / m_vortex_radius) *
                floe::geometry::direct_orthogonal(vortex_center_to_pt);
        } else {
            return m_vortex_max_norm * (m_vortex_radius / distance_to_center) *
                floe::geometry::direct_orthogonal(vortex_center_to_pt) / distance_to_center;
        }
    }

};

template <typename TPoint>
TPoint
PhysicalData<TPoint>::water_speed(point_type pt) {
    point_type resp;
    if (m_water_mode == 0){
        resp = geostrophic_water_speed();
    } else {
        resp = get_speed(pt, m_water_mode);
    }
    return resp + m_geo_relative_water_speed;
}

template <typename TPoint>
TPoint
PhysicalData<TPoint>::air_speed(point_type pt) {
    if (m_air_mode == 0){
        return topaz_air_speed(pt);
    } else {
        return get_speed(pt, m_air_mode);
    }
}

template <typename TPoint>
void
PhysicalData<TPoint>::load_matlab_topaz_data(std::string const& filename){
    floe::io::matlab::read_topaz_from_file(filename, m_ocean_data_hours, m_air_data_hours);
    interpolate_hour_to_minute();
}


template <typename TPoint>
void
PhysicalData<TPoint>::interpolate_hour_to_minute()
{
    interpolate_hour_to_minute(m_ocean_data_hours, m_ocean_data_minutes);
    interpolate_hour_to_minute(m_air_data_hours, m_air_data_minutes);
}


template <typename TPoint>
void
PhysicalData<TPoint>::interpolate_hour_to_minute(point_vector const& data_hours, point_vector& data_minutes){
    if (!data_hours.size())
        return;

    point_type p0, p1, P0, P1, pm, Pm;

    // init phase
    const real_type init_hours = 12;
    p0 = {0,0};
    P0 = {norm2(p0), atan2(p0.y, p0.x)}; // polar coordinates of P0.
    p1 = data_hours[0];
    P1 = {norm2(p1), atan2(p1.y, p1.x)};
    for (std::size_t j = 0; j!= 60 * init_hours; ++j)
    {
        real_type h_frac = (real_type)j/ (60 * init_hours);
        if (std::abs(P1[1] - P0[1]) > M_PI) // to avoid doing more than a U turn
            P0[1] += copysign(2 * M_PI, P1[1] - P0[1]);
        Pm = (1. - h_frac) * P0 + h_frac * P1;
        pm = {Pm[0] * cos(Pm[1]), Pm[0] * sin(Pm[1])};
        data_minutes.push_back(pm);
    }

    // we store polar coordinates in cartesian point, to get simple interpolation on norm and angle.
    P0 = P1; // polar coordinates of P0.
    for (std::size_t i = 1; i != data_hours.size(); ++i)
    {
        p1 = data_hours[i];
        P1 = {norm2(p1), atan2(p1.y, p1.x)};
        for (std::size_t j = 0; j!=60; ++j)
        {
            real_type h_frac = (real_type)j/60;
            if (std::abs(P1[1] - P0[1]) > M_PI) // to avoid doing more than a U turn
                P0[1] += copysign(2 * M_PI, P1[1] - P0[1]);
            Pm = (1. - h_frac) * P0 + h_frac * P1;
            pm = {Pm[0] * cos(Pm[1]), Pm[0] * sin(Pm[1])};
            data_minutes.push_back(pm);
        }
        P0 = P1;
    }
}


template <typename TPoint>
TPoint
PhysicalData<TPoint>::geostrophic_water_speed(point_type p){
    return minute_value(m_time_ref, m_ocean_data_minutes);
}


template <typename TPoint>
TPoint
PhysicalData<TPoint>::topaz_air_speed(point_type p){
    return minute_value(m_time_ref, m_air_data_minutes);
}


template <typename TPoint>
TPoint
PhysicalData<TPoint>::minute_value(real_type t, point_vector const& data_minutes){
    std::size_t minutes = t / 60;
    if (minutes < data_minutes.size())
        return data_minutes[minutes];
    else
        return data_minutes[-1];
}


template <typename TPoint>
void
PhysicalData<TPoint>::update_water_speed(point_type diff_speed)
{
    m_geo_relative_water_speed += diff_speed;
}

template <typename TPoint>
void
PhysicalData<TPoint>::init_random_vortex(){
    real_type avg_Rc{300000}, delta_Rc{110000};
    real_type avg_V{8}, delta_V{4};
    real_type min_Umax{15}, max_Umax{30};
    auto dist_Rc = std::uniform_real_distribution<double>{
        avg_Rc - delta_Rc, avg_Rc + delta_Rc};
    auto dist_V = std::uniform_real_distribution<double>{
        avg_V - delta_V, avg_V + delta_V};
    auto dist_Umax = std::uniform_real_distribution<double>{
        min_Umax, max_Umax};
    auto dist_angle = std::uniform_real_distribution<double>{
        0, 2 * M_PI};
    auto gen = std::default_random_engine{};
    std::cout << "RANDOM SEED " << std::time(0);
    // gen.seed(std::time(0));
    gen.seed(0);
    m_vortex_radius = dist_Rc(gen);
    m_vortex_max_norm = dist_Umax(gen);
    real_type theta = dist_angle(gen);
    m_vortex_origin = 1e6 * point_type{cos(theta), sin(theta)}; // Vortex starts at 1000km from origin
    m_vortex_speed = - dist_V(gen) * point_type{cos(theta), sin(theta)};
    std::cout << "WIND VORTEX : "
    << "radius = " << m_vortex_radius
    << ", speed =  " << m_vortex_speed
    << ", origin = " << m_vortex_origin
    << ", max norm = " << m_vortex_max_norm << std::endl;
}

template <typename TPoint>
TPoint
PhysicalData<TPoint>::get_speed(point_type pt, int mode){
    point_type resp;
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


}} // namespace floe::dynamics


#endif // OPE_PHYSICAL_DATA_HPP