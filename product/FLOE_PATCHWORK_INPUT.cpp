#include <iostream>
#include <cassert>
#include "../product/config/config.hpp"
#include <string>
#include <vector>
#include <random>
#include <ctime>
// Boost geometry
#include "floe/geometry/frame/frame_transformers.hpp"

/*
Create an unit hdf5 input file (domain = 1m * 1m) from a list of input files
by randomly placing floe packs on a N * P grid
*/

template<typename T>
using scale_transformer = boost::geometry::strategy::transform::scale_transformer<T, 2, 2>;

int main( int argc, char* argv[] )
{   
    using namespace std;
    using namespace types;

    if ( argc < 4 )
    {
        cout << "Usage: " << argv[0] << " <nb_rows> <nb_cols> <input_file_name_1> <input_file_name_2> ...<input_file_name_N>" << endl;
        return 1;
    }

    int nb_input_files = argc - 3;
    std::cout << nb_input_files << " files" << std::endl;
    int nb_rows = atoi(argv[1]);
    int nb_cols = atoi(argv[2]);
    std::vector<std::string> input_filenames;
    for (int i = 0; i < nb_input_files; i++)
        input_filenames.push_back(argv[3+i]);
    for (auto& f : input_filenames)
        std::cout << f << std::endl;

    std::default_random_engine generator;
    generator.seed(std::time(0)); // seed for not having pseudorandom
    std::uniform_int_distribution<int> rand_int_distrib(0, nb_input_files - 1);
    std::vector<int> floe_packs_limits;
    problem_type P;
    auto& floes = P.get_floe_group().get_floes();
    int id_floe = 0;
    // decltype(w) nxp_unit_window{-nb_cols/2, nb_cols/2, -nb_rows/2, nb_rows/2};
    for (int i=0; i<nb_rows; i++){
        for (int j=0; j<nb_cols; j++){
            point_type trans{
                (value_type)(-nb_cols + 2*j + 1) / 2,
                (value_type)(-nb_rows + 2*i + 1) / 2
            };
            int idx = rand_int_distrib(generator);
            std::string& input_filename = input_filenames[idx];
            P.load_config(input_filename);
            auto w = P.get_floe_group().get_initial_window();
            point_type win_center{(w[0]+w[1])/2, (w[2]+w[3])/2};
            auto win_width = (w[1] - w[0]);
            floe_packs_limits.push_back(floes.size());
            for (id_floe; id_floe < floes.size(); id_floe++){
                auto & floe = floes[id_floe];
                // re-center pack
                floe.state().pos -= win_center;
                // Resize input (to 1x1 unit domain)
                auto base_shape = floe.static_floe().geometry();
                floe::geometry::transform( base_shape, floe.static_floe().geometry(), scale_transformer<value_type>{ 1/win_width } );
                point_type base_pos = floe.state().pos;
                floe::geometry::transform( base_pos, floe.state().pos, scale_transformer<value_type>{ 1/win_width } );
                // translate floe
                floe.state().pos += trans;
            }
        }
    }

    P.get_floe_group().set_initial_window({(value_type)-nb_cols/2, (value_type)nb_cols/2, (value_type)-nb_rows/2, (value_type)nb_rows/2});
    P.make_input_file();
    return 0;
}