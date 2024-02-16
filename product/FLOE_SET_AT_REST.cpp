#include <iostream>
#include <cassert>
#include "../product/config/interrupt.hpp"
#include "../product/config/config.hpp"
#include "../product/config/config_base.hpp"


/*
Loads an input file, stops floes in the defined window, and writes a new input file.
*/


int main( int argc, char* argv[] )
{   
    using namespace std;
    using namespace types;

    if ( argc < 2 )
    {
        cout << "Usage: " << argv[0] << " <input_file_name> " << endl;
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
    // bool generate_floes = false;
    std::string matlab_topaz_filename = "io/library/DataTopaz01.mat";

    problem_type P;
    P.QUIT = &QUIT;
    P.load_config(input_list_floe_filename);

    std::cout << "Kinetic energy before setting at rest: " << P.get_floe_group().kinetic_energy() << "\n";

    std::array<value_type, 4>    window_type    = P.get_floe_group().get_initial_window();
    value_type                   win_width      = 2*(window_type[1]-window_type[0]);
    value_type                   win_height     = 2*(window_type[3]-window_type[2]);
    std::cout << "size of the window: " << win_width/2 << " | " << win_height/2 << "\n";
    P.get_floe_group().stop_floes_in_window(win_width, win_height);

    std::cout << "Kinetic energy after setting at rest: " << P.get_floe_group().kinetic_energy() << std::endl;
    
    P.make_input_file();

    return 0;
}