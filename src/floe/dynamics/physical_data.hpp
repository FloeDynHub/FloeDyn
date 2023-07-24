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
#include <cassert>


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
        m_firstVortexZoneDistToOrigin{0}, m_vortexZoneSize{0}, m_nbVortexByZone{0}, m_nb_vortex{0},
        m_vortex_radius{}, m_vortex_origin{}, m_vortex_speed{}, m_vortex_max_norm{}, m_nb_time_step{}, m_dt{300},
        m_water_mode{0}, m_air_mode{0} {}

    //! water speed accessor (m/s)
    point_type water_speed(point_type pt = {0,0});
    //! air speed accessor (m/s)
    point_type air_speed(point_type pt = {0,0});
    //! OBL update speed (m/s)
    void update_water_speed(point_type diff_speed);
    //! OBL speed accessor for output (m/s)
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
    void set_modes(int air_mode, int water_mode) {
        m_air_mode      = air_mode;
        m_water_mode    = water_mode;

        if (m_air_mode==0 && m_water_mode==0) {std::cout << "NO Atmospheric and Ocean currents!" << std::endl;}
        else if (m_air_mode==1 && m_water_mode==1) {std::cout << "Atmospheric and Ocean currents from weather data files" << std::endl;}
        else if ((m_air_mode==2 && m_water_mode==0) || (m_air_mode==0 && m_water_mode==2)) {
            std::cout << "Convergent current towards the target area for initial ice pack generation" << std::endl;
        }
        else if ((m_air_mode==4 && m_water_mode==0) || (m_air_mode==0 && m_water_mode==4)) {
            std::cout << "Current directed towards the obstacle" << std::endl;
        }
        else if (m_air_mode==5) {
            this->init_vortex();m_water_mode = 0;
            std::cout << "Storm defined as a wind vortex" << std::endl;
        }
        else if (m_air_mode==6) {
            this->init_random_vortex();m_water_mode = 0;
            std::cout << "Storm defined as a wind vortex" << std::endl;
        }
        else { std::cout << "Error: air and/or water modes: " << m_air_mode << " and " << m_water_mode << " are unknown!" << std::endl; assert(true==false); }
    }
    
    //!< Air and water speeds:
    void set_speeds(real_type air_speed, real_type water_speed) {
        m_air_speed     = air_speed;
        m_water_speed   = water_speed;

        if ((m_air_mode==2 && m_water_mode==0) || (m_air_mode==0 && m_water_mode==2) 
            || (m_air_mode==4 && m_water_mode==0) || (m_air_mode==0 && m_water_mode==4)) {
            std::cout << "The Atmospheric current speed is fixed to: " << m_air_speed << "m/s   |   The Ocean current speed is fixed to: "
                << m_water_speed << "m/s" << std::endl;
        }
    }
    //!< Initialization of the storm as a wind vortex
    void init_random_vortex();
    void init_vortex();

    //!< For checking the vortex wind evolution
    real_type get_vortex_wind_speed(std::size_t i) {return set_vortex_wind_speed(i);};

    //!< air mode accessor 
    int get_air_mode() {return m_air_mode;};


    //!< vortex getter and setter
    size_t get_nb_vortex() const {return m_nb_vortex;};
    void set_nb_vortex(size_t nb_vortex) {  m_nb_vortex = nb_vortex;  };
    size_t get_nbVortexByZone() const {return m_nbVortexByZone;};
    void set_nbVortexByZone(size_t nbVortexByZone) {m_nbVortexByZone = nbVortexByZone;};
    real_type get_vortexZoneSize() const {return m_vortexZoneSize;};
    void set_vortexZoneSize(real_type vortexZoneSize) {m_vortexZoneSize = vortexZoneSize;};
    real_type get_firstVortexZoneDistToOrigin() const {return m_firstVortexZoneDistToOrigin;};
    void set_firstVortexZoneDistToOrigin(real_type firstVortexZoneDistToOrigin) {m_firstVortexZoneDistToOrigin = firstVortexZoneDistToOrigin;};

private:

    point_vector m_ocean_data_hours; //!< Geostrophic datas
    point_vector m_air_data_hours; //!< Geostrophic datas
    point_vector m_ocean_data_minutes; //!< Geostrophic datas
    point_vector m_air_data_minutes; //!< Geostrophic datas
    real_type const& m_time_ref; //!< reference to time variable in second

    point_type m_geo_relative_water_speed; //!< Water speed correction compared to geostrophic data

    // window dimension (for generator)
    real_type m_window_width;
    real_type m_window_height;
    // vortex attributes
    real_type                m_firstVortexZoneDistToOrigin; // in [m]
    real_type                m_vortexZoneSize; // in [m]
    std::size_t              m_nbVortexByZone; 
    std::size_t              m_nb_vortex;
    std::vector<real_type>   m_vortex_radius;
    point_vector             m_vortex_origin;    //!< vortex eye
    point_vector             m_vortex_speed;
    std::vector<real_type>   m_vortex_max_norm;
    std::vector<int>         m_nb_time_step; // number of velocity changement before reaching Umax for a vortex
    real_type   m_dt;         //!< time in second for the velocity discretization.
    // modes
    int m_water_mode;
    int m_air_mode;
    //!< speeds
    real_type m_air_speed;
    real_type m_water_speed;

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
    point_type get_speed(point_type pt = {0,0}, int mode=0, real_type speed=0);
    //! center convergent current (negative coeff will give divergent field)
    // point_type centered_convergent_field(point_type pt = {0,0}, real_type coeff = 1) {
    //     return - coeff * pt / norm2(pt);
    // }
    //! convergent current outside null rectangle window
    point_type convergent_outside_window_field(point_type pt = {0,0}, real_type speed=1) {
        real_type x{0}, y{0};
        if (std::abs(pt.x) > m_window_width / 2) x = - speed * pt.x / std::abs(pt.x);
        if (std::abs(pt.y) > m_window_height / 2) y = - speed * pt.y / std::abs(pt.y);
        return {x,y};
    }
    //! circular current inside the windows field
    // point_type circular_inside_window_field(point_type pt = {0,0}, real_type speed_current=1) {
    //     real_type x{0}, y{0};
    //     if (std::abs(pt.x) <= m_window_width / 2 and std::abs(pt.y) <= m_window_height / 2){
    //         if (std::abs(pt.x) < std::abs(pt.y)) x = - 10 * pt.y / std::abs(pt.y);
    //         if (std::abs(pt.x) >= std::abs(pt.y)) y = 10 * pt.x / std::abs(pt.x);
    //     }
    //     return {x,y};
    // }
    // ! convergent current outside null rectangle window
    // point_type convergent_outside_window_circular_inside_field(
    //     point_type pt = {0,0}, real_type out_speed=30, real_type in_speed=10
    // ){
    //     return convergent_outside_window_field(pt, out_speed)
    //         + circular_inside_window_field(pt, in_speed);
    // }
    //! x axis convergent wind, then constant

    point_type x_convergent_then_constant(point_type pt = {0,0}, real_type coeff=1, real_type speed=1)
    {
        real_type x{0}, y{0};
        if (m_time_ref < 500) {
            if (pt.x<=-0.3) {x=1e-3; y=-5e-2;}
        }
        else {x=speed;}
        return {x,y};        
        // if (m_time_ref < 1500) {
        //     return {0, - coeff * (pt.y / (100 + std::abs(pt.y)))};
        // }
        // else {
        //     point_type speed_current = {speed,0};
        //     return speed_current;
        // }
    }
    //! vortex storm
    point_type vortex_center(std::size_t i){
        return m_vortex_origin[i] + m_time_ref * m_vortex_speed[i];
    }

    //!< set the vortex wind speed:
    real_type set_vortex_wind_speed(std::size_t i) {
        real_type vortex_wind_speed;

        if ( m_time_ref < (m_nb_time_step[i]+1)*m_dt ) {
            vortex_wind_speed = int(m_time_ref/m_dt) * m_vortex_max_norm[i]/m_nb_time_step[i];
        }
        else if ( m_time_ref < 2*m_nb_time_step[i]*m_dt ) {
            vortex_wind_speed = int( ( 2*m_nb_time_step[i]*m_dt - m_time_ref )/m_dt ) *
                m_vortex_max_norm[i]/m_nb_time_step[i];
        }
        else {vortex_wind_speed=0;}

        return vortex_wind_speed;
    }

    point_type vortex(point_type pt) {
        point_type totalAirVelocityAppliedToPoint{0,0};

        bool display=false;

        for (std::size_t i=0; i<m_nb_vortex; ++i) {
            point_type vortex_center_to_pt = pt - vortex_center(i);
            real_type distance_to_center = norm2(vortex_center_to_pt);
            point_type unitOrthoDirection = floe::geometry::direct_orthogonal(vortex_center_to_pt)/distance_to_center;

            //!< wind velocity computation:
            real_type vortex_wind_speed = set_vortex_wind_speed(i);

            if (distance_to_center < m_vortex_radius[i]){
                totalAirVelocityAppliedToPoint += vortex_wind_speed * (distance_to_center / m_vortex_radius[i]) *
                    unitOrthoDirection;
            } else {
                totalAirVelocityAppliedToPoint += vortex_wind_speed * std::pow(m_vortex_radius[i] / distance_to_center,4) *
                    unitOrthoDirection;
            }

            // if (vortex_wind_speed>0) 
            // {
            //     display = true;
            //     std::cout << "WIND VORTEX number " << i
            //     << ", distToCenter = " << distance_to_center
            //     << ", orthoDir =  " << unitOrthoDirection
            //     << ", speed = " << vortex_wind_speed
            //     << ", radius = " << m_vortex_radius[i]
            //     << ", totalVel = " << totalAirVelocityAppliedToPoint << std::endl;
            // }

        }


        // if (display) {std::cout << "point: " << pt << std::endl;}
        return totalAirVelocityAppliedToPoint;
    }

};

template <typename TPoint>
TPoint
PhysicalData<TPoint>::water_speed(point_type pt) {
    point_type resp;
    if (m_water_mode == 1){
        resp = geostrophic_water_speed();
    } else {
        resp = get_speed(pt, m_water_mode, m_water_speed);
    }
    return resp + m_geo_relative_water_speed;
}

template <typename TPoint>
TPoint
PhysicalData<TPoint>::air_speed(point_type pt) {
    if (m_air_mode == 1){
        return topaz_air_speed(pt);
    } else {
        return get_speed(pt, m_air_mode, m_air_speed);
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
PhysicalData<TPoint>::init_vortex(){

	m_nb_vortex=3;
	m_vortex_radius={291720,277476,342577};//291720; //277476; //342577; //325639; //374848; //dist_Rc(gen);
	m_vortex_max_norm={23.562,21.4638,27.6211};//23.562; //21.4638; //27.6211; //19.1004; //19.3233; //dist_Umax(gen);
	m_vortex_origin={point_type{-529615 , 848238},point_type{932306 , -361672},point_type{-377083 , -926180}};// Vortex starts at 1000km from origin
	m_vortex_speed={point_type{3.00029 , -4.8053},point_type{-8.7237 , 3.38421},point_type{3.73322 , 9.16942}}; //point_type{1.26598 , 8.87528}; //point_type{-9.10191 , 0.813075}; //- dist_V(gen) * point_type{cos(theta), sin(theta)};

	
    std::cout << "how many vortex?: " << m_nb_vortex << std::endl;

    if (m_nb_vortex==0) {m_nb_time_step.push_back(1);}

    for (std::size_t i=0; i<m_nb_vortex; ++i) {  

        //!< linear increase of the vortex wind speed from 0 to m_vortex_max_norm reached above the ice pack:
        real_type dist_to_ice_pack_center = norm2(m_vortex_origin[i]);
        m_nb_time_step.push_back( int( dist_to_ice_pack_center / (norm2(m_vortex_speed[i]) * m_dt) ) );
        std::cout << "m_dn " << m_nb_time_step[i] << std::endl;
        if (m_nb_time_step[i]<=0) {std::cout << "Problem with the distance to the ice pack and the velocity of the eye vortex." << std::endl;}
        assert(m_nb_time_step[i]>0);

        std::cout << "WIND VORTEX : "
        << "radius = " << m_vortex_radius[i]
        << ", velocity =  " << m_vortex_speed[i]
        << ", speed = " << norm2(m_vortex_speed[i])
        << ", origin = " << m_vortex_origin[i]
        << ", max norm = " << m_vortex_max_norm[i]
        << ", the vortex takes: " << int(m_nb_time_step[i]*m_dt/3600) << "h to reach the ice pack center.\n" << std::endl;
    }
}


template <typename TPoint>
void
PhysicalData<TPoint>::init_random_vortex(){

    real_type avg_Rc{5000}, delta_Rc{3000};
    real_type delta_dist{20000};
    real_type avg_V{8}, delta_V{4};
    real_type min_Umax{15}, max_Umax{30};
    auto dist_Rc = std::uniform_real_distribution<double>{
        avg_Rc - delta_Rc, avg_Rc + delta_Rc};
    auto dist_delta = std::uniform_real_distribution<double>{
        0, delta_dist};
    auto dist_V = std::uniform_real_distribution<double>{
        avg_V - delta_V, avg_V + delta_V};
    auto dist_Umax = std::uniform_real_distribution<double>{
        min_Umax, max_Umax};
    auto dist_angle = std::uniform_real_distribution<double>{
        0, 2 * M_PI};
    auto dist_quartAngle = std::uniform_real_distribution<double>{
        M_PI/2, M_PI};
    auto delay = std::uniform_real_distribution<double>{
        -40000, 40000};
    auto gen = std::default_random_engine{};

    std::size_t k=0;

    // creer boucle pour definir une serie de vortex!! On peut definir les m_vortex_... charac comme
    // des std::vector<real_type> ou des points_vector=std::vector>point_type> qui contient
    // les centre des vortex, leur vitesse.
    // dans la function vortex on peut boucler sur les vortex et sommer les vents pour le calcul des vents.

    //!< warning    with MPI simulation, std::time(0) may be different from workers!!
    // std::cout << "RANDOM SEED " << std::time(0);
    // gen.seed(90839527656);
     
    gen.seed(645214440);
    real_type theta_base = dist_angle(gen);
    real_type theta_total = theta_base;
    std::cout << "angle de base: " << theta_base << std::endl;
    // gen.seed(0);

    std::cout << "how many vortex?: " << m_nb_vortex << std::endl;
    std::cout << "how many vortex by zone?: " << m_nbVortexByZone << std::endl;

    if (m_nb_vortex==0) {m_nb_time_step.push_back(1);}

    for (std::size_t i=0; i<m_nb_vortex; ++i) {
        m_vortex_radius.push_back(dist_Rc(gen)); //291720; //277476; //342577; //325639; //374848; //dist_Rc(gen);
        m_vortex_max_norm.push_back(dist_Umax(gen)); //23.562; //21.4638; //27.6211; //19.1004; //19.3233; //dist_Umax(gen);
        if (i>0) {theta_total += dist_quartAngle(gen);}

        if (i>=(k+1)*m_nbVortexByZone) {k+=1;}
        real_type distToOrigin = m_firstVortexZoneDistToOrigin + k*m_vortexZoneSize + dist_delta(gen);
        std::cout << "distToOrigin: " << distToOrigin << std::endl;

        m_vortex_origin.push_back(distToOrigin * point_type{cos(theta_total), sin(theta_total)}); //point_type{-529615 , 848238}; //point_type{932306 , -361672}; //point_type{-377083 , -926180}; //point_type{-141212 , -989979}; //point_type{996034 , -88975.9}; //1e6 * point_type{cos(theta), sin(theta)}; // Vortex starts at 1000km from origin
        std::cout << "angle supplementaire: " << theta_total << " | cos(theta total): " << cos(theta_total) << std::endl;

        real_type delayX = delay(gen);
        real_type delayY = delay(gen);
        std::cout << "delay: (" << delayX << "," << delayY << ")" << std::endl;
        point_type vectorToDelay = point_type{delayX-m_vortex_origin[i].x , delayY-m_vortex_origin[i].y};
        point_type velDirToOrigin_withDelay = dist_V(gen) * vectorToDelay / norm2(vectorToDelay);
        std::cout << "vectorToDelay: " << vectorToDelay << std::endl;
        m_vortex_speed.push_back(velDirToOrigin_withDelay); //point_type{3.00029 , -4.8053}; //point_type{-8.7237 , 3.38421}; //point_type{3.73322 , 9.16942}; //point_type{1.26598 , 8.87528}; //point_type{-9.10191 , 0.813075}; //- dist_V(gen) * point_type{cos(theta), sin(theta)};

        //!< linear increase of the vortex wind speed from 0 to m_vortex_max_norm reached above the ice pack:
        real_type dist_to_ice_pack_center = norm2(m_vortex_origin[i]);
        m_nb_time_step.push_back( int( dist_to_ice_pack_center / (norm2(m_vortex_speed[i]) * m_dt) ) );
        if (m_nb_time_step[i]<=0) {std::cout << "Problem with the distance to the ice pack and the velocity of the eye vortex." << std::endl;}
        assert(m_nb_time_step[i]>0);

        std::cout << "WIND VORTEX : "
        << "radius = " << m_vortex_radius[i]
        << ", velocity =  " << m_vortex_speed[i]
        << ", speed = " << norm2(m_vortex_speed[i])
        << ", origin = " << m_vortex_origin[i]
        << ", max norm = " << m_vortex_max_norm[i]
        << ", the vortex takes: " << int(m_nb_time_step[i]*m_dt/3600) << "h to reach the ice pack center.\n" << std::endl;
    }
}



template <typename TPoint>
TPoint
PhysicalData<TPoint>::get_speed(point_type pt, int mode, real_type speed){
    point_type resp;
    switch(mode){
        // case 6:
        //     // std::cout << "generation of a centered convergent field\n";
        //     return centered_convergent_field(pt, speed);
        case 2:
            // std::cout << "generation of a convergent outside window field\n";
            return convergent_outside_window_field(pt, speed);
        // case 3:
        //     // std::cout << "generation of a convergent outside window circular inside field\n";
        //     return convergent_outside_window_circular_inside_field(pt, speed);
        case 4:
            // std::cout << "generation of a x convergent then constant\n";
            return x_convergent_then_constant(pt, speed);
        case 5:
            // std::cout << "generation of a vortex storm\n";
            return vortex(pt);
        case 0:
            return {0,0};

        default :  
            return {0,0};   

    }
}

}} // namespace floe::dynamics


#endif // OPE_PHYSICAL_DATA_HPP
