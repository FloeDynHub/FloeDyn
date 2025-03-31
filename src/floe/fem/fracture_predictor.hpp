/**
 * @file fracture_predictor.hpp
 * @author Silouane
 * @brief stores features of a floe and its impacts to predict fracture using a specified model
 * @version 0.1
 * @date 2025-01-15
 *
 *
 */

#ifndef FRACTURE_PREDICTOR
#define FRACTURE_PREDICTOR

#include <onnxruntime_cxx_api.h>
#include <cpu_provider_factory.h>

namespace floe { namespace fem {

/**
 * @brief fast prediction of fracture. 
 * 
 * @details This class is used to predict the fracture of a floe using a pre-trained model. It computes the corresponding features and predicts the probability of fracture using a pre-trained model. The features may be send out to the log in order to build a training database.
 * 
 * @tparam TFloe 
 */
template <
    typename TFloe
>
class FracturePredictor
{
public :
    using real_type = typename TFloe::real_type;
    using point_type = typename TFloe::point_type;
    typedef TFloe floe_type;
    typedef typename floe_type::mesh_type     mesh_type;
    using multi_point_type = typename mesh_type::multi_point_type;
    using connectivity_type = typename mesh_type::connectivity_type;

    FracturePredictor();
    bool prepare_entry_floe(floe_type * floe);
    bool prepare_entry_impact(size_t i_impact, point_type impact);
    bool clear_entry_impact();
    std::string get_database_entry() {return m_entry_floe.str() + m_entry_impact.str();};
    std::vector<double> caliper_diameters(std::vector<point_type> geometry, size_t n_angles);
    std::vector<double> calculate_moments(std::vector<point_type> geometry);
    bool predict_fracture();
    bool prepare_predictor();

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

    Ort::Env m_env;
    Ort::SessionOptions m_session_options;
    std::unique_ptr<Ort::Session> m_session;
    std::vector<std::string> m_input_node_names;
    std::vector<std::string> m_output_node_names;

    float m_probaThreshold;
    size_t m_batch_size;
    size_t m_num_classes;
    std::vector<int64_t> m_input_dims;
    std::vector<int64_t> m_output_dims; 
    std::string m_model_path; 

};


template <typename TFloe>
FracturePredictor<TFloe>::FracturePredictor():
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
    m_entry_impact{std::stringstream("")},
    m_env(ORT_LOGGING_LEVEL_WARNING, "onnx_example"),
    m_session_options(),
    m_probaThreshold{0.1},
    m_batch_size{1},
    m_num_classes{2},
    m_input_dims{std::vector<int64_t>({1, 19})}, // Batch size = 1, 19 features
    m_output_dims{std::vector<int64_t>({1, 2})}, // 1 sample, 2 possible classes T/F.
    m_model_path("./fracture_predictors/fracture_predictor_RF.onnx")
{
    prepare_predictor();
}

template <typename TFloe>
std::vector<double> FracturePredictor<TFloe>::caliper_diameters(std::vector<point_type> geometry, size_t n_angles)
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
std::vector<double> FracturePredictor<TFloe>::calculate_moments(std::vector<point_type> geometry)
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
bool FracturePredictor<TFloe>::prepare_entry_floe(floe_type * floe)
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
bool FracturePredictor<TFloe>::clear_entry_impact()
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
bool FracturePredictor<TFloe>::prepare_entry_impact(size_t i_impact, point_type impact)
{
    m_entry_impact.str("");
    m_entry_impact.clear();
    m_impact_x += impact.x;
    m_impact_y += impact.y;
    m_impact += norm2(impact);
    // là, éventuellement, faudra plus tard ajouter une rotation pour se mettre dans le repère du floe. Non en fait pas sûr ! Ptêt ben que non en fait ! à vérifier avec le setup one ball one beam mais ch'pense que pas b'soin. 
    m_phi_0 += m_mincd_x * impact.x + m_mincd_y * impact.y;
    m_phi_1 += m_maxcd_x * impact.x + m_maxcd_y * impact.y;
    m_phi_2 += (m_mincd_x * impact.x + m_mincd_y * impact.y) / m_mincd;
    m_phi_3 += (m_maxcd_x * impact.x + m_maxcd_y * impact.y) / m_maxcd;
    m_phi_4 += (m_mincd_x * impact.x + m_mincd_y * impact.y) / (m_mincd * m_mincd);
    m_phi_5 += (m_maxcd_x * impact.x + m_maxcd_y * impact.y) / (m_maxcd * m_maxcd);
    m_entry_impact << "; impact_id = " << i_impact << " ; impact_y = " << m_impact_y << " ; impact_x = " << m_impact_x << " ; impact = " << m_impact << " ; phi_0 = " << m_phi_0 << " ; phi_1 = " << m_phi_1 << " ; phi_2 = " << m_phi_2 << " ; phi_3 = " << m_phi_3 << " ; phi_4 = " << m_phi_4 << " ; phi_5 = " << m_phi_5 << " ; ";
    return true;
}

template <typename TFloe>
bool FracturePredictor<TFloe>::prepare_predictor()
{
    // std::string model_path = "/Users/silouane/Documents/code/Floe/FloeDyn/bazar_silou/piml/model_RF.onnx";
    // std::string model_path = "./fracture_predictor/model_RF.onnx";

    m_session_options.SetIntraOpNumThreads(1);
    m_session_options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
    
    try {
        m_session = std::make_unique<Ort::Session>(m_env, m_model_path.c_str(), m_session_options);
    } catch (const Ort::Exception& e) {
        std::cerr << "Failed to load ONNX model: " << e.what() << std::endl;
        return false;
    }

    Ort::AllocatorWithDefaultOptions allocator;
    size_t num_inputs = m_session->GetInputCount();
    size_t num_outputs = m_session->GetOutputCount();

    std::vector<std::string> input_node_names_temp;
    std::vector<std::string> output_node_names_temp;

    input_node_names_temp.resize(num_inputs);
    output_node_names_temp.resize(num_outputs);
    m_input_node_names.resize(num_inputs);
    m_output_node_names.resize(num_outputs);

    // storing names in strings because otherwise they are destroyed at the end of the loop
    for (size_t i = 0; i < num_inputs; ++i) {
        auto input_name = m_session->GetInputNameAllocated(i, allocator);
        input_node_names_temp[i] = input_name.get();
        m_input_node_names[i] = input_node_names_temp[i].c_str();
    }

    for (size_t i = 0; i < num_outputs; ++i) {
        auto output_name = m_session->GetOutputNameAllocated(i, allocator);
        output_node_names_temp[i] = output_name.get();
        m_output_node_names[i] = output_node_names_temp[i].c_str();
    }

    return true;
}


template <typename TFloe>
bool FracturePredictor<TFloe>::predict_fracture()
{
    std::cout << "Predicting fracture with PIML" << std::endl;

    if (!m_session) {
        std::cout << "ONNX model not correctly loaded. Make sure you call prepare_predictor() first." << std::endl;
        return true;
    }

    Ort::MemoryInfo memory_info = Ort::MemoryInfo::CreateCpu(
        OrtArenaAllocator, OrtMemTypeDefault);

    // Inference on this current sample
    std::vector<double> input_data_double = {m_area, m_perimeter, m_surf_over_length, m_surf_over_length_squared, m_impact_y, m_impact_x, m_I_x, m_I_y, m_J_x, m_J_y, m_mincd, m_maxcd, m_meancd, m_phi_0, m_phi_1, m_phi_2, m_phi_3, m_phi_4, m_phi_5};
    std::vector<float> input_data = std::vector<float>(input_data_double.begin(), input_data_double.end());
    
    Ort::Value input_tensor = Ort::Value::CreateTensor<float>(
        memory_info, input_data.data(), input_data.size(),
        m_input_dims.data(), m_input_dims.size());
    
    std::vector<float> output_probabilities(1 * m_num_classes);
    try {
        //The Ort::Value stores the output
        Ort::Value output_tensor = Ort::Value::CreateTensor<float>(
            memory_info,
            output_probabilities.data(),
            output_probabilities.size(),
            m_output_dims.data(),
            m_output_dims.size());

        // Inference
        try {
            const char* output_probability_name = m_output_node_names[1].c_str(); // backtransforming from string to const char* 
            std::vector<const char*> input_node_names = {m_input_node_names[0].c_str()}; // idem 
            m_session->Run(
                Ort::RunOptions{nullptr},
                input_node_names.data(), &input_tensor, 1, // inouts
                &output_probability_name, &output_tensor, 1   // outputs
            );
        } catch (const Ort::Exception& e) {
            std::cerr << "Error during inference: " << e.what() << std::endl;
            return 1;
        }
        // for (size_t i = 0; i < m_batch_size; ++i) 
        // {
        //     std::cout << "Probabilities (T F): ";
        //     for (size_t j = 0; j < m_num_classes; ++j) {
        //         std::cout << output_probabilities[i * m_num_classes + j] << " ";
        //     }
        //     std::cout << std::endl;
        // }
    } catch (const Ort::Exception& e) {
        std::cerr << "Error during inference: " << e.what() << std::endl;
        return 1;
    }

    if (output_probabilities[0] > m_probaThreshold)
    {
        std::cout << "The floe might break (p = " << output_probabilities[0] << " )" << std::endl;
        return true;
    }
    else
    {
        std::cout << "The floe should not break (p = " << output_probabilities[0] << " )" << std::endl;
        // std::cout << "Predictor says it should not break" << std::endl;
        return false; 
    }


    return true;
}


}} // namespace floe::fem

#endif
