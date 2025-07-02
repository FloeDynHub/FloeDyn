#include <iostream>
#include <cassert>
#include "../product/config/config.hpp"
#include <limits>
#include <math.h>
#include <boost/version.hpp>

/*
Reads an input file and gives some informations on the floe pack
*/

int main( int argc, char* argv[] )
{   
    using namespace std;
    using namespace types;

    std::cout << "Using Boost "     
          << BOOST_VERSION / 100000     << "."  // major version
          << BOOST_VERSION / 100 % 1000 << "."  // minor version
          << BOOST_VERSION % 100                // patch level
          << std::endl;

    if ( argc < 2 )
    {
        cout << "Usage: " << argv[0] << " <input_file_name>" << endl;
        return 1;
    }

    std::string input_list_floe_filename = argv[1];
    problem_type P;
    P.load_config(input_list_floe_filename);
    auto& floes = P.get_floe_group().get_floes();
    auto w = P.get_floe_group().get_initial_window();
    auto win_width = (w[1] - w[0]);
    auto win_height = (w[3] - w[2]);
    std::cout << "nb_floes :: " << P.get_floe_group().nb_floes() << std::endl;
    std::cout << "concentration :: " << (int)(P.get_floe_group().initial_concentration() * 100) << std::endl;
    std::cout << "Domain_w_h :: " << win_width << "," << win_height << std::endl;
    std::cout << "initial_window :: " << w[0] << "," << w[1] << "," << w[2] << "," << w[3] << std::endl;
    auto wc = P.get_floe_group().bounding_window(0);
    std::cout << "current_window :: " << wc[0] << "," << wc[1] << "," << wc[2] << "," << wc[3] << std::endl;
    std::cout << "current_concentration :: " << P.get_floe_group().floe_concentration() << std::endl;
    std::cout << "kinetic_energy :: " << P.get_floe_group().kinetic_energy() << std::endl;

    value_type max_radius = 0;
    value_type min_radius = std::numeric_limits<value_type>::max();
    P.create_optim_vars();
    for (auto optim_ptr : P.proximity_detector().data().get_optims()){
        value_type floe_global_disk_radius = optim_ptr->global_disk().radius;
        max_radius = std::max(max_radius, floe_global_disk_radius);
        min_radius = std::min(min_radius, floe_global_disk_radius);
    }
    std::cout << "floe_diameter_max :: " << max_radius * 2 << std::endl;
    std::cout << "floe_diameter_min :: " << min_radius * 2 << std::endl;
    // MPI grid dim constraint : max_floe_radius < cell_size / 4;
    int n = 1;
    while (win_width / (n+1) > 4 * max_radius) n++; 
    std::cout << "mpi_max_grid_dim :: " << n << std::endl;
    std::cout << "mpi_average_floes_per_core " << floes.size() / (n*n) << std::endl;

    return 0;
}