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
    sigaction(SIGTERM,&sa,NULL);
    sigaction(SIGSEGV,&sa,NULL);
    sigaction(SIGINT,&sa,NULL); 
    sigaction(SIGILL,&sa,NULL);
    sigaction(SIGABRT,&sa,NULL);
    sigaction(SIGFPE,&sa,NULL); 


    DT_DEFAULT = atof(argv[3]);

    std::string matlab_list_floe_filename = argv[1];
    std::string matlab_topaz_filename = "/Users/jouetq/Projects/FloeCpp/io/DataTopaz01.mat";

    problem_type P;
    P.load_matlab_config(matlab_list_floe_filename);
    P.load_matlab_topaz_data(matlab_topaz_filename);
    P.auto_topology();

    if (argc == 6)
        P.recover_states_from_file("/Users/jouetq/Projects/FloeCpp/io/in.h5", atof(argv[5]));
    

    P.solve(atoi(argv[2]), atof(argv[4]));

    return 0;
}