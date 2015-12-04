double DT_DEFAULT;
#include <iostream>
#include <cassert>
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
    // omp_set_num_threads(1);
    Eigen::initParallel();
    #endif

    DT_DEFAULT = atof(argv[3]);

    std::string matlab_list_floe_filename = argv[1];
    bool generate_floes = false;
    if (matlab_list_floe_filename == "generator") generate_floes = true; 
    std::string matlab_topaz_filename = "io/inputs/DataTopaz01.mat";

    problem_type P;
    if (!generate_floes)
        P.load_matlab_config(matlab_list_floe_filename);
    else {
        generator_type G;
        std::cout << "How many floes ? ";
        int nb_floes;
        std::cin >> nb_floes;
        std::cout << "Which concentration (between 0 and 1) ? ";
        double concentration;
        std::cin >> concentration;
        G.generate_floe_set(nb_floes, concentration);
        P.set_floe_group(G.get_floe_group());
        #ifdef PBC
        auto win = G.get_window();
        P.set_topology(topology_type(win[0], win[1], win[2], win[3]));
        #endif
    }

    std::cout << "read TOPAZ" << std::endl;
    P.load_matlab_topaz_data(matlab_topaz_filename);


    if (argc == 7)
    {
        std::cout << "Please enter an input filename : ";
        std::string strFilename;
        std::cin >> strFilename;
        P.recover_states_from_file("io/outputs/" + strFilename, atof(argv[6]));
    }

    
    std::cout << "SOLVE..." << std::endl;
    P.solve(atoi(argv[2]), DT_DEFAULT, atoi(argv[5]), atof(argv[4]));

    return 0;
}