/*!
 * \file ope/physical_data.hpp
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


namespace floe { namespace ope
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
    using value_type = decltype(TPoint::x);
    using point_vector = std::vector<point_type>;

    //! Constructor
    PhysicalData(value_type const& time_ref) :
        m_ocean_data_hours{}, m_air_data_hours{},
        m_ocean_data_minutes{}, m_air_data_minutes{},
        m_time_ref{time_ref}, m_water_speed{0,0},
        m_geo_relative_water_speed{0, 0} {}

    //! water speed accessor
    point_type water_speed(point_type = {0,0});
    //! air speed accessor
    point_type air_speed(point_type = {0,0});
    void update_water_speed(point_type diff_speed);
    //! OBL speed accessor for output
    inline point_type OBL_speed() const { return m_geo_relative_water_speed; }

    //! Load ocean and wind data from a topaz file
    void load_matlab_topaz_data(std::string const& filename);

private:

    point_vector m_ocean_data_hours; //!< Geostrophic datas
    point_vector m_air_data_hours; //!< Geostrophic datas
    point_vector m_ocean_data_minutes; //!< Geostrophic datas
    point_vector m_air_data_minutes; //!< Geostrophic datas
    value_type const& m_time_ref; //!< reference to time variable

    point_type m_water_speed; //!< Corrected water speed
    point_type m_geo_relative_water_speed; //!< Water speed correction compared to geostrophic data

    //! Interpolate all hour datas to minute datas
    void interpolate_hour_to_minute();
    //! Interpolate individual hour datas to minute datas
    void interpolate_hour_to_minute(point_vector const&, point_vector&);
    //! Get value in minute datas from time
    point_type minute_value(value_type t, point_vector const&);
    //! Get geostrophic water speed
    point_type geostrophic_water_speed(point_type = {0,0});

};


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
    const value_type init_hours = 12;
    p0 = {0,0};
    P0 = {norm2(p0), atan2(p0.y, p0.x)}; // polar coordinates of P0.
    p1 = data_hours[0];
    P1 = {norm2(p1), atan2(p1.y, p1.x)};
    for (std::size_t j = 0; j!= 60 * init_hours; ++j)
    {
        value_type h_frac = (value_type)j/ (60 * init_hours);
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
            value_type h_frac = (value_type)j/60;
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
PhysicalData<TPoint>::water_speed(point_type p){
    return m_water_speed;
}


template <typename TPoint>
TPoint
PhysicalData<TPoint>::air_speed(point_type p){
    return minute_value(m_time_ref, m_air_data_minutes);
}


template <typename TPoint>
TPoint
PhysicalData<TPoint>::minute_value(value_type t, point_vector const& data_minutes){
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
    m_water_speed = geostrophic_water_speed() + m_geo_relative_water_speed;
}



}} // namespace floe::ope


#endif // OPE_PHYSICAL_DATA_HPP