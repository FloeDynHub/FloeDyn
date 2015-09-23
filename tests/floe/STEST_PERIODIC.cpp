#include <iostream>
#include "../tests/floe/interrupt.hpp"
#include "../tests/floe/config_periodic.hpp"


int main( int argc, char* argv[] )
{   
    using namespace std;

    if ( argc < 5 )
    {
        cout << "Usage: " << argv[0] << " <matlab_file_name> <end_time> <dt_default> <out_step_nb>" << endl;
        return 1;
    }

    // Handling interruption signal
    struct sigaction sa;
    memset( &sa, 0, sizeof(sa) );
    sa.sa_handler = interruption::got_signal;
    sigfillset(&sa.sa_mask);
    sigaction(SIGINT,&sa,NULL); 
    // sigaction(SIGABRT,&sa,NULL);
    // sigaction(SIGTERM,&sa,NULL);

    // sigaction(SIGSEGV,&sa,NULL);
    // sigaction(SIGILL,&sa,NULL);
    // sigaction(SIGFPE,&sa,NULL); 

    #ifdef _OPENMP
    // omp_set_num_threads(0);
    Eigen::initParallel();
    #endif

    DT_DEFAULT = atof(argv[3]);

    std::string matlab_list_floe_filename = argv[1];
    std::string matlab_topaz_filename = "io/DataTopaz01.mat";

    problem_type P;
    P.load_matlab_config(matlab_list_floe_filename);
    P.load_matlab_topaz_data(matlab_topaz_filename);
    // P.auto_topology();
    P.set_topology(floe::io::matlab::read_pze_from_file<topology_type>(matlab_list_floe_filename));

    if (argc == 6)
        P.recover_states_from_file("io/in.h5", atof(argv[5]));
    

    P.solve(atoi(argv[2]), atof(argv[4]));

    return 0;
}