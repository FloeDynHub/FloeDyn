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

// #include <onnxruntime_cxx_api.h>
// #include <cpu_provider_factory.h>

#include "floe/fem/fracture_descriptor.hpp"

namespace floe { namespace fem {

/**
 * @brief fast prediction of fracture. 
 * 
 * @details This class is used to infere the fracture of a floe using a pre-trained model. The features come from a FractureDescriptor instance. The model is loaded in the prepare() stage, call in the constructor of the PartialFloeGroup. The class uses ONNX runtime to perform the inference.
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
    bool predict_fracture(fem::FractureDescriptor<floe_type> const & descriptor) const;
    bool prepare_predictor();

private :
    // Ort::Env m_env;
    // Ort::SessionOptions m_session_options;
    // std::unique_ptr<Ort::Session> m_session;
    // std::vector<std::string> m_input_node_names;
    // std::vector<std::string> m_output_node_names;

    // float m_probaThreshold;
    // size_t m_batch_size;
    // size_t m_num_classes;
    // std::vector<int64_t> m_input_dims;
    // std::vector<int64_t> m_output_dims; 
    // std::string m_model_path; 

};


template <typename TFloe>
FracturePredictor<TFloe>::FracturePredictor()//:
    // m_env(ORT_LOGGING_LEVEL_WARNING, "onnx_example"),
    // m_session_options(),
    // m_probaThreshold{0.1},
    // m_batch_size{1},
    // m_num_classes{2},
    // m_input_dims{std::vector<int64_t>({1, 19})}, // Batch size = 1, 19 features
    // m_output_dims{std::vector<int64_t>({1, 2})}, // 1 sample, 2 possible classes (T/F).
    // m_model_path("./fracture_predictors/fracture_predictor_RF.onnx")
{}

template <typename TFloe>
bool FracturePredictor<TFloe>::prepare_predictor()
{
    // m_session_options.SetIntraOpNumThreads(1);
    // m_session_options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
    // try {
    //     m_session = std::make_unique<Ort::Session>(m_env, m_model_path.c_str(), m_session_options);
    // } catch (const Ort::Exception& e) {
    //     std::cerr << "Failed to load ONNX model: " << e.what() << std::endl;
    //     return false;
    // }

    // Ort::AllocatorWithDefaultOptions allocator;
    // size_t num_inputs = m_session->GetInputCount();
    // size_t num_outputs = m_session->GetOutputCount();

    // std::vector<std::string> input_node_names_temp;
    // std::vector<std::string> output_node_names_temp;

    // input_node_names_temp.resize(num_inputs);
    // output_node_names_temp.resize(num_outputs);
    // m_input_node_names.resize(num_inputs);
    // m_output_node_names.resize(num_outputs);

    // // storing names in strings because otherwise they are destroyed at the end of the loop
    // for (size_t i = 0; i < num_inputs; ++i) {
    //     auto input_name = m_session->GetInputNameAllocated(i, allocator);
    //     input_node_names_temp[i] = input_name.get();
    //     m_input_node_names[i] = input_node_names_temp[i].c_str();
    // }

    // for (size_t i = 0; i < num_outputs; ++i) {
    //     auto output_name = m_session->GetOutputNameAllocated(i, allocator);
    //     output_node_names_temp[i] = output_name.get();
    //     m_output_node_names[i] = output_node_names_temp[i].c_str();
    // }

    return true;
}


template <typename TFloe>
bool FracturePredictor<TFloe>::predict_fracture(fem::FractureDescriptor<floe_type> const & descriptor) const
{
    // std::cout << "Predicting fracture with PIML" << std::endl;

    // if (!m_session) {
    //     std::cout << "ONNX model not correctly loaded. Make sure prepare_predictor() is called first." << std::endl;
    //     return true;
    // }

    // Ort::MemoryInfo memory_info = Ort::MemoryInfo::CreateCpu(
    //     OrtArenaAllocator, OrtMemTypeDefault);

    // // Inference on this current sample
    // std::vector<double> input_data_double = {
    //     descriptor.area(), 
    //     descriptor.perimeter(), 
    //     descriptor.surf_over_length(), 
    //     descriptor.surf_over_length_squared(),
    //     descriptor.impact_y(), 
    //     descriptor.impact_x(), 
    //     descriptor.I_x(),
    //     descriptor.I_y(),
    //     descriptor.J_x(),
    //     descriptor.J_y(), 
    //     descriptor.mincd(), 
    //     descriptor.maxcd(), 
    //     descriptor.meancd(), 
    //     descriptor.phi_0(),
    //     descriptor.phi_1(),
    //     descriptor.phi_2(), 
    //     descriptor.phi_3(), 
    //     descriptor.phi_4(),
    //     descriptor.phi_5()};
    
    // // std::vector<double> input_data_double = {1.83463443e+04,  4.97598508e+02,  3.68697735e+01, 7.40954262e-02,  7.11362000e-04,  6.42058000e-05, 1.19617137e+08,  1.04247884e+08,  2.99902673e+07, 2.52056500e+07,  1.43067902e+02,  1.70361736e+02, 1.58263600e+02,  3.36889734e-02,  1.21568467e-01, 2.35475413e-04,  7.13590209e-04,  1.64589967e-06, 4.18867654e-06};

    // std::vector<float> input_data = std::vector<float>(input_data_double.begin(), input_data_double.end());
    
    // Ort::Value input_tensor = Ort::Value::CreateTensor<float>(
    //     memory_info, input_data.data(), input_data.size(),
    //     m_input_dims.data(), m_input_dims.size());
    
    // std::vector<float> output_probabilities(1 * m_num_classes);
    // try {
    //     //The Ort::Value stores the output
    //     Ort::Value output_tensor = Ort::Value::CreateTensor<float>(
    //         memory_info,
    //         output_probabilities.data(),
    //         output_probabilities.size(),
    //         m_output_dims.data(),
    //         m_output_dims.size());

    //     // Inference
    //     try {
    //         const char* output_probability_name = m_output_node_names[1].c_str(); // backtransforming from string to const char* 
    //         std::vector<const char*> input_node_names = {m_input_node_names[0].c_str()}; // idem 
    //         m_session->Run(
    //             Ort::RunOptions{nullptr},
    //             input_node_names.data(), &input_tensor, 1, // inouts
    //             &output_probability_name, &output_tensor, 1   // outputs
    //         );
    //     } catch (const Ort::Exception& e) {
    //         std::cerr << "Error during inference: " << e.what() << std::endl;
    //         return 1;
    //     }
    // } catch (const Ort::Exception& e) {
    //     std::cerr << "Error during inference: " << e.what() << std::endl;
    //     return 1;
    // }

    // if (output_probabilities[0] > m_probaThreshold)
    // {
    //     std::cout << "The floe might break (p = " << output_probabilities[0] << " )" << std::endl;
    //     return true;
    // }
    // else
    // {
    //     std::cout << "The floe should not break (p = " << output_probabilities[0] << " )" << std::endl;
    //     return false; 
    // }
    return true;
}


}} // namespace floe::fem

#endif
