#include <iostream>
#include <cassert>
#include "../product/config/interrupt.hpp"
#include "../product/config/config.hpp"
#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;

/* Function used to check that of 'for_what' is specified, then
   'required_option' is specified too. */
void option_dependency(const po::variables_map& vm,
                        const char* for_what, const char* required_option)
{
    if (vm.count(for_what) && !vm[for_what].defaulted())
        if (vm.count(required_option) == 0 || vm[required_option].defaulted())
            throw logic_error(std::string("Option '") + for_what 
                              + "' requires option '" + required_option + "'.");
}

void handle_exception(std::exception& e){
    std::cerr << "Error: " << e.what() << "\n";
}

int main( int argc, char* argv[] )
{   
    using namespace types;

    string input_file_name;
    value_type endtime;
    value_type default_time_step = 10;
    value_type out_time_step = 60;
    int OBL_status = 0;
    value_type epsilon = 0.4;
    string matlab_topaz_filename = "io/inputs/DataTopaz01.mat";
    
    po::options_description desc("Allowed options");
    desc.add_options()
    // First parameter describes option name/short name
    // The second is parameter to option
    // The third is description
    ("help,h", "print usage message")
    ("input,i", po::value(&input_file_name)->required(), "input file path")
    ("tend,t", po::value(&endtime)->required(), "simulation duration (seconds)")
    ("step,s", po::value(&default_time_step), "default time step")
    ("outstep,o", po::value(&out_time_step), "output time step")
    ("obl", po::value(&OBL_status), "OBL status (0 or 1)")
    ("rectime,r", po::value<value_type>(), "time to recover states from")
    ("recfile,f", po::value<string>(), "file name to recover states from")
    ("nbfloes,n", po::value<int>(), "generator : how many floes ?")
    ("concentration,c", po::value<value_type>(), "generator : floes concentration (between 0 and 1)")
    ("epsilon,e", po::value(&epsilon), "collision restitution coeff")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {  
        cout << desc << "\n";
        return 0;
    }

    try {
        po::notify(vm);
        option_dependency(vm, "rectime", "recfile");
      } catch (std::exception& e) {
        handle_exception(e);
        return 1;
      }

    // Handling interruption signal
    struct sigaction sa; 
    memset( &sa, 0, sizeof(sa) );
    sa.sa_handler = interruption::got_signal;
    sigfillset(&sa.sa_mask);
    sigaction(SIGINT,&sa,NULL);
    sigaction(SIGSEGV,&sa,NULL); 

    #ifdef _OPENMP
    // omp_set_num_threads(1);
    Eigen::initParallel();
    #endif

    bool generate_floes = false;
    if (input_file_name == "generator") generate_floes = true; 

    problem_type P(epsilon, OBL_status);
    P.QUIT = &QUIT;
    if (!generate_floes){
        try {
            P.load_config(input_file_name);
        }
        catch(std::exception& e)
        {
            handle_exception(e);
            return 1;
        }
    }
    
    for (auto& floe : P.get_floe_group().floes()){
        auto& shape = floe.static_floe().geometry();
        auto& mesh = floe.static_floe().mesh();
        auto mass_center = floe::integration::integrate(
            [] (value_type x, value_type y) { return point_type{x, y}; },
            mesh, integration_strategy()
        ) / floe::integration::integrate(
            [] (value_type x, value_type y) { return 1.; },
            mesh,integration_strategy()
        );
        // std::cout << "MC : " << mass_center << std::endl;
        polygon_type shape_cpy = shape;
        geometry::transform( shape_cpy, shape, geometry::frame::transformer( typename floe_type::frame_type{-mass_center, 0} ));
        mesh_type mesh_cpy = mesh;
        geometry::transform( mesh_cpy, mesh, geometry::frame::transformer( typename floe_type::frame_type{-mass_center, 0} ));
        floe.state().x += mass_center.x;
        floe.state().y += mass_center.y;
    }
    // cout.precision(17);
    // std::cout << P.get_floe_group().total_area();

    return 0;
}