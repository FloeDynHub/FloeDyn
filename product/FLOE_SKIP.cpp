#include <iostream>
#include <cassert>
#include "../product/config/interrupt.hpp"
#include "../product/config/config.hpp"


int main( int argc, char* argv[] )
{   
    using namespace std; 

    if ( argc < 7 )
    {
        cout << "Usage: " << argv[0] << " <matlab_file_name> <end_time> <dt_default> <out_step_nb> <OBL status> <nbskip>" << endl;
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

    FLOE_SKIP = atoi(argv[6]);
    std::cout << "#SKIP : " << FLOE_SKIP << std::endl;

    problem_type P;
    P.load_matlab_config(matlab_list_floe_filename);

    std::cout << "read TOPAZ" << std::endl;
    P.load_matlab_topaz_data(matlab_topaz_filename);

    std::cout << "SOLVE..." << std::endl;
    P.solve(atoi(argv[2]), atof(argv[3]);, atoi(argv[5]), atof(argv[4]));


        
    return 0;
}