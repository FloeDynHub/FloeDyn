#include <iostream>
#include "../product/interrupt.hpp"
#include "../product/config.hpp"


int main( int argc, char* argv[] )
{   
    using namespace std;

    if ( argc < 6 )
    {
        cout << "Usage: " << argv[0] << " <matlab_file_name> <end_time> <dt_default> <out_step_nb> <OBL status>" << endl;
        return 1;
    }

    // Handling interruption signal
    struct sigaction sa;
    memset( &sa, 0, sizeof(sa) );
    sa.sa_handler = interruption::got_signal;
    sigfillset(&sa.sa_mask);
    sigaction(SIGINT,&sa,NULL); 

    #ifdef _OPENMP
    // omp_set_num_threads(0);
    Eigen::initParallel();
    #endif

    DT_DEFAULT = atof(argv[3]);

    std::string matlab_list_floe_filename = argv[1];
    std::string matlab_topaz_filename = "io/DataTopaz01.mat";

    problem_type P;
    P.load_matlab_config(matlab_list_floe_filename);
    std::cout << "read TOPAZ" << std::endl;
    P.load_matlab_topaz_data(matlab_topaz_filename);

    OBL_STATUS = atoi(argv[5]);

    if (argc == 7)
    {
        std::cout << "Please enter an input filename : ";
        std::string strFilename;
        std::cin >> strFilename;
        P.recover_states_from_file("io/" + strFilename, atof(argv[6]));
    }
    
    std::cout << "SOLVE..." << std::endl;
    P.solve(atoi(argv[2]), atof(argv[4]));

    return 0;
}