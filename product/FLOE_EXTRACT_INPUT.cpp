#include <iostream>
#include <cassert>
#include "../product/config/interrupt.hpp"
#include "../product/config/config.hpp"

/*
Create an hdf5 input file from output file
*/


int main( int argc, char* argv[] )
{   
    using namespace std;
    using namespace types;

    if ( argc < 2 )
    {
        cout << "Usage: " << argv[0] << " <input_file_name> <rec_time>" << endl;
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

    std::string input_list_floe_filename = argv[1];
    bool generate_floes = false;
    std::string matlab_topaz_filename = "io/library/DataTopaz01.mat";

    problem_type P;
    P.QUIT = &QUIT;
    P.load_config(input_list_floe_filename);
    if (argc > 2){
        double rectime = atof(argv[2]);
        P.load_matlab_topaz_data(matlab_topaz_filename); // segfault otherwise (?)
        P.recover_states_from_file(input_list_floe_filename, rectime, false);
    }
    P.make_input_file();
    // cout << P.get_floe_group().total_area();

    return 0;
}