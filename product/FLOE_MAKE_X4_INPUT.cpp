#include <iostream>
#include <cassert>
#include "../product/config/config.hpp"

/*
Create an hdf5 input file from an input file by loading it 4 times and translating floes.
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

    std::string input_list_floe_filename = argv[1];

    problem_type P;
    auto& floes = P.get_floe_group().get_floes();
    P.load_config(input_list_floe_filename);
    int nb_floes = floes.size();
    auto w = P.get_floe_group().get_initial_window();
    auto X = (w[1] - w[0]) / 2;
    auto Y = (w[3] - w[2]) / 2;
    std::vector<types::point_type> trans_vector;
    trans_vector.push_back({-X, -Y});
    P.load_config(input_list_floe_filename);
    trans_vector.push_back({-X, Y});
    P.load_config(input_list_floe_filename);
    trans_vector.push_back({X, Y});
    P.load_config(input_list_floe_filename);
    trans_vector.push_back({X, -Y});
    for (int i = 0; i < 4; ++i){
        for (int id = i * nb_floes; id < (i+1) * nb_floes; ++id){
            auto& floe = floes[id];
            floe.state().pos += trans_vector[i];
        }
    }
    P.get_floe_group().set_initial_window({{ w[0] - X, w[1] + X, w[2] - Y, w[3] + Y }});
    P.make_input_file();

    return 0;
}