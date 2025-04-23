/**
 * @file fracture_descriptor.hpp
 * @author Silouane
 * @brief stores features of a floe and its impacts to predict fracture using a specified model
 * @version 0.1
 * @date 2025-01-15
 *
 *
 */

#ifndef FRACTURE_DESCRIPTOR
#define FRACTURE_DESCRIPTOR


namespace floe { namespace fem {

/**
 * @brief Feature storage for fracture prediction. 
 * 
 * @details This class stores data about the floe and its last impact. It contains all features used as input for the fast predictors. The features may be send out to the log in order to build a training database, or given to the predictor to predict the fracture of the floe.
 * 
 * @tparam TFloe 
 */
template <
    typename TFloe
>
class FractureDescriptor
{
public :
    using real_type = typename TFloe::real_type;
    using point_type = typename TFloe::point_type;
    typedef TFloe floe_type;
    typedef typename floe_type::mesh_type     mesh_type;
    using multi_point_type = typename mesh_type::multi_point_type;
    using connectivity_type = typename mesh_type::connectivity_type;

    FractureDescriptor();
    bool prepare_entry_floe(floe_type * floe);
    bool prepare_entry_impact(size_t i_impact, point_type impact);
    bool clear_entry_impact();
    std::string get_database_entry() {return m_entry_floe.str() + m_entry_impact.str();};
    std::vector<double> caliper_diameters(std::vector<point_type> geometry, size_t n_angles);
    std::vector<double> calculate_moments(std::vector<point_type> geometry);

    inline double area() const {return m_area;}
    inline double perimeter() const {return m_perimeter;}
    inline double surf_over_length() const {return m_surf_over_length;}
    inline double surf_over_length_squared() const {return m_surf_over_length_squared;}
    inline double I_x() const {return m_I_x;}
    inline double I_y() const {return m_I_y;}
    inline double J_x() const {return m_J_x;}
    inline double J_y() const {return m_J_y;}
    inline double mincd() const {return m_mincd;}
    inline double mincd_x() const {return m_mincd_x;}
    inline double mincd_y() const {return m_mincd_y;}
    inline double maxcd() const {return m_maxcd;}
    inline double maxcd_x() const {return m_maxcd_x;}
    inline double maxcd_y() const {return m_maxcd_y;}
    inline double meancd() const {return m_meancd;}
    inline double impact_y() const {return m_impact_y;}
    inline double impact_x() const {return m_impact_x;}
    inline double impact() const {return m_impact;}
    inline double phi_0() const {return m_phi_0;}
    inline double phi_1() const {return m_phi_1;}
    inline double phi_2() const {return m_phi_2;}
    inline double phi_3() const {return m_phi_3;}
    inline double phi_4() const {return m_phi_4;}
    inline double phi_5() const {return m_phi_5;}

private :
    double m_area;
    double m_perimeter;
    double m_surf_over_length;
    double m_surf_over_length_squared;
    double m_I_x;
    double m_I_y;
    double m_J_x;
    double m_J_y;
    double m_mincd;
    double m_mincd_x;
    double m_mincd_y;
    double m_maxcd;
    double m_maxcd_x;
    double m_maxcd_y;
    double m_meancd;
    double m_impact_y;
    double m_impact_x;
    double m_impact;
    double m_phi_0;
    double m_phi_1;
    double m_phi_2;
    double m_phi_3;
    double m_phi_4;
    double m_phi_5;
    std::stringstream m_entry_floe;
    std::stringstream m_entry_impact;
};


template <typename TFloe>
FractureDescriptor<TFloe>::FractureDescriptor():
    m_area{0},
    m_perimeter{0},
    m_surf_over_length{0},
    m_surf_over_length_squared{0},
    m_I_x{0},
    m_I_y{0},
    m_J_x{0},
    m_J_y{0},
    m_mincd{0},
    m_mincd_x{0},
    m_mincd_y{0},
    m_maxcd{0},
    m_maxcd_x{0},
    m_maxcd_y{0},
    m_meancd{0},
    m_impact_y{0},
    m_impact_x{0},
    m_impact{0},
    m_phi_0{0},
    m_phi_1{0},
    m_phi_2{0},
    m_phi_3{0},
    m_phi_4{0},
    m_phi_5{0},
    m_entry_floe{std::stringstream("")},
    m_entry_impact{std::stringstream("")}
{}

template <typename TFloe>
std::vector<double> FractureDescriptor<TFloe>::caliper_diameters(std::vector<point_type> geometry, size_t n_angles)
{
    double mincd_x = 0;
    double mincd_y = 0;
    double maxcd_x = 0;
    double maxcd_y = 0;
    double mincd = 1e18;
    double maxcd = 0;
    double diameter_sum = 0.0;

    std::vector<double> caliper_diameters;
    caliper_diameters.reserve(7);
    
    std::vector<point_type> points = geometry;
    size_t n_points = points.size();
    if (n_points < 2)
    {
        std::cerr << "Geometry must have at least two points." << std::endl;
        return caliper_diameters;
    }

    for (size_t i = 0; i < n_angles; ++i)
    {
        double angle = i * M_PI / n_angles;
        double cos_angle = std::cos(angle);
        double sin_angle = std::sin(angle);
        double max_x(0);
        double min_x(1e18);
        for (size_t j = 0; j < n_points; ++j)
        {
            double x = points[j].x;
            double y = points[j].y;
            double x_rot = x * cos_angle - y * sin_angle;
            if (x_rot > max_x)
                max_x = x_rot;
            if (x_rot < min_x)
                min_x = x_rot;
        }
        double diameter = max_x - min_x;
        diameter_sum += diameter;
        if (diameter > maxcd)
        {
            maxcd = diameter;
            maxcd_x = -diameter * cos_angle;
            maxcd_y = diameter * sin_angle;
        }
        if (diameter < mincd)
        {
            mincd = diameter;
            mincd_x = -diameter * cos_angle;
            mincd_y = diameter * sin_angle;
        }
    }
    double meancd = diameter_sum / n_angles;
    caliper_diameters.push_back(mincd);
    caliper_diameters.push_back(mincd_x);
    caliper_diameters.push_back(mincd_y);
    caliper_diameters.push_back(maxcd);
    caliper_diameters.push_back(maxcd_x);
    caliper_diameters.push_back(maxcd_y);
    caliper_diameters.push_back(meancd);
    return caliper_diameters;
}


template <typename TFloe>
std::vector<double> FractureDescriptor<TFloe>::calculate_moments(std::vector<point_type> geometry)
{
    double J_x = 0;
    double J_y = 0;
    double J_xy = 0;
    double I_x = 0;
    double I_y = 0;
    double I_xy = 0;

    std::vector<double> moments;
    moments.reserve(6);

    std::vector<point_type> points = geometry;
    size_t n_points = points.size();
    if (n_points < 2)
    {
        std::cerr << "Geometry must have at least two points." << std::endl;
        return moments;
    }

    double centroid_x = 0;
    double centroid_y = 0;
    for (size_t i = 0; i < n_points; ++i)
    {
        centroid_x += points[i].x;
        centroid_y += points[i].y;
    }
    centroid_x /= n_points;
    centroid_y /= n_points;

    for (size_t i = 0; i < n_points; ++i)
    {
        double x0 = points[i].x - centroid_x;
        double y0 = points[i].y - centroid_y;
        double x1 = points[(i + 1) % n_points].x - centroid_x;
        double y1 = points[(i + 1) % n_points].y - centroid_y;

        double a = x0 * y1 - x1 * y0;

        J_x += a * (y0 * y0 + y0 * y1 + y1 * y1) / 12;
        J_y += a * (x0 * x0 + x0 * x1 + x1 * x1) / 12;
        J_xy += a * (x0 * y1 + 2 * x0 * y0 + 2 * x1 * y1 + x1 * y0) / 24;

        I_x += a * y0 * y0;
        I_y += a * x0 * x0;
        I_xy += a * x0 * y0;
    }

    moments.push_back(J_x);
    moments.push_back(J_y);
    moments.push_back(J_xy);
    moments.push_back(I_x);
    moments.push_back(I_y);
    moments.push_back(I_xy);
    return moments;
}


template <typename TFloe>
bool FractureDescriptor<TFloe>::prepare_entry_floe(floe_type * floe)
{
    m_entry_floe.str("");
    m_entry_floe.clear();
    
    auto geometry = floe->static_floe().geometry();
    std::vector<point_type> boundary = geometry.outer();

    m_area = floe->static_floe().area();
    m_perimeter = boost::geometry::perimeter(geometry);
    m_surf_over_length = m_area / m_perimeter;
    m_surf_over_length_squared = m_surf_over_length * m_surf_over_length;
    std::vector<double> moments = calculate_moments(boundary);
    std::vector<double> caliper_diam = caliper_diameters(boundary, 180);
    m_I_x = moments[0];
    m_I_y = moments[1];
    m_J_x = moments[3];
    m_J_y = moments[4];
    m_mincd = caliper_diam[0];
    m_mincd_x = caliper_diam[1];
    m_mincd_y = caliper_diam[2];
    m_maxcd =  caliper_diam[3]; 
    m_maxcd_x = caliper_diam[4];
    m_maxcd_y = caliper_diam[5];
    m_meancd = caliper_diam[6];
    m_entry_floe << "area = " << m_area << " ; perimeter = " << m_perimeter << " ; surf_over_length = " << m_surf_over_length << " ; surf_over_length_squared = " << m_surf_over_length_squared << " ; I_x = " << m_I_x << " ; I_y = " << m_I_y << " ; J_x = " << m_J_x << " ; J_y = " << m_J_y << " ; mincd = " << m_mincd << " ; mincd_x = " << m_mincd_x << " ; mincd_y = " << m_mincd_y << " ; maxcd = " << m_maxcd << " ; maxcd_x = " << m_maxcd_x << " ; maxcd_y = " << m_maxcd_y << " ; meancd = " << m_meancd << " ; ";
    return true;
}

template <typename TFloe>
bool FractureDescriptor<TFloe>::clear_entry_impact()
{
    m_entry_impact.str("");
    m_entry_impact.clear();
    m_impact_x = 0;
    m_impact_y = 0;
    m_impact = 0;
    m_phi_0 = 0;
    m_phi_1 = 0;
    m_phi_2 = 0;
    m_phi_3 = 0;
    m_phi_4 = 0;
    m_phi_5 = 0;
    return true;
}

template <typename TFloe>
bool FractureDescriptor<TFloe>::prepare_entry_impact(size_t i_impact, point_type impact)
{
    m_entry_impact.str("");
    m_entry_impact.clear();
    m_impact_x += impact.x;
    m_impact_y += impact.y;
    m_impact += norm2(impact);
    // bug potentiel : faut-il ajouter une rotation pour se mettre dans le repère du floe ? à vérifier avec le setup one ball one beam 
    m_phi_0 += m_mincd_x * impact.x + m_mincd_y * impact.y;
    m_phi_1 += m_maxcd_x * impact.x + m_maxcd_y * impact.y;
    m_phi_2 += (m_mincd_x * impact.x + m_mincd_y * impact.y) / m_mincd;
    m_phi_3 += (m_maxcd_x * impact.x + m_maxcd_y * impact.y) / m_maxcd;
    m_phi_4 += (m_mincd_x * impact.x + m_mincd_y * impact.y) / (m_mincd * m_mincd);
    m_phi_5 += (m_maxcd_x * impact.x + m_maxcd_y * impact.y) / (m_maxcd * m_maxcd);
    m_entry_impact << "; impact_id = " << i_impact << " ; impact_y = " << m_impact_y << " ; impact_x = " << m_impact_x << " ; impact = " << m_impact << " ; phi_0 = " << m_phi_0 << " ; phi_1 = " << m_phi_1 << " ; phi_2 = " << m_phi_2 << " ; phi_3 = " << m_phi_3 << " ; phi_4 = " << m_phi_4 << " ; phi_5 = " << m_phi_5 << " ; ";
    return true;
}



}} // namespace floe::fem

#endif
