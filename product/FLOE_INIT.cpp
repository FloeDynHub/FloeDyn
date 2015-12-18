#include <iostream>
#include <cassert>
#include "../product/config/interrupt.hpp"
#include "../product/config/config.hpp"


int main( int argc, char* argv[] )
{   
    using namespace std; 

    if ( argc < 7 )
    {
        cout << "Usage: " << argv[0] << " <matlab_file_name> <end_time> <dt_default> <out_step_nb> <OBL status> <time init>" << endl;
        return 1;
    }

    // Handling interruption signal
    struct sigaction sa; 
    memset( &sa, 0, sizeof(sa) );
    sa.sa_handler = interruption::got_signal;
    sigfillset(&sa.sa_mask);
    sigaction(SIGINT,&sa,NULL); 

    #ifdef _OPENMP
    // omp_set_num_threads(1);
    Eigen::initParallel();
    #endif

    std::string matlab_list_floe_filename = argv[1];
    bool generate_floes = false;
    if (matlab_list_floe_filename == "generator") generate_floes = true; 
    std::string matlab_topaz_filename = "io/inputs/DataTopaz01.mat";

    types::generator_problem_type P1;
    P1.load_matlab_config(matlab_list_floe_filename);
    auto& init_physical_data = P1.get_dynamics_manager().get_external_forces().get_physical_data();
    init_physical_data.set_window_size(0,0);
    init_physical_data.set_modes(2,0);
    std::cout << "INIT..." << std::endl;
    P1.solve(atoi(argv[2]), atof(argv[3]), atoi(argv[5]), atof(argv[4]));

    types::problem_type P;
    std::cout << "read TOPAZ" << std::endl;
    P.load_matlab_topaz_data(matlab_topaz_filename);
    P.set_floe_group(P1.get_floe_group());
    std::cout << "SOLVE..." << std::endl;
    P.solve(atoi(argv[2]), atof(argv[3]), atoi(argv[5]), atof(argv[4]));
        
    return 0;
}