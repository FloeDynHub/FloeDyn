#include <iostream>
#include <cassert>
#include "../product/config/config.hpp"

// Boost geometry
#include "floe/geometry/frame/frame_transformers.hpp"
#include <cmath>
#include <cstddef>     // std::size_t
#include <type_traits> // std::enable_if
#include <ostream>
#include "floe/arithmetic/container_operators.hpp"

/*
Create an hdf5 input file from an input file by applying an homothety to floes' geometries.
*/

template<typename T>
using scale_transformer = boost::geometry::strategy::transform::scale_transformer<T, 2, 2>;

int main( int argc, char* argv[] )
{   
    using namespace std;
    using namespace types;

    if ( argc < 3 )
    {
        cout << "Usage: " << argv[0] << " <input_file_name> <window_width (m)>" << endl;
        return 1;
    }

    std::string input_list_floe_filename = argv[1];
    int new_win_width = atof(argv[2]);

    problem_type P;
    auto& floes = P.get_floe_group().get_floes();
    P.load_config(input_list_floe_filename);
    auto w = P.get_floe_group().get_initial_window();
    auto win_width = (w[1] - w[0]);
    // auto win_height = (w[3] - w[2]);
    value_type scale_coeff = new_win_width / win_width;
    decltype(w) new_window{
        {scale_coeff * w[0],
        scale_coeff * w[1],
        scale_coeff * w[2],
        scale_coeff * w[3]}
    };
    for (auto& floe : floes){
        auto base_shape = floe.static_floe().geometry();
        floe::geometry::transform( base_shape, floe.static_floe().geometry(), scale_transformer<value_type>{ scale_coeff } );
        point_type base_pos = floe.state().pos;
        floe::geometry::transform( base_pos, floe.state().pos, scale_transformer<value_type>{ scale_coeff } );
    }
    P.get_floe_group().set_initial_window(new_window);
    P.make_input_file();

    return 0;
}