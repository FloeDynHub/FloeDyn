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
        cout << "Usage: " << argv[0] << " <input_file_name> <x trans (m)> <y trans (m)>" << endl;
        return 1;
    }

    std::string input_list_floe_filename = argv[1];
    types::point_type trans = {atof(argv[2]), atof(argv[3])};

    problem_type P;
    auto& floes = P.get_floe_group().get_floes();
    P.load_config(input_list_floe_filename);
    auto w = P.get_floe_group().get_initial_window();
    decltype(w) new_window{
        {w[0] + trans[0],
         w[1] + trans[1],
         w[2] + trans[0],
         w[3] + trans[1]}    
    };
    for (auto& floe : floes){
        floe.state().pos += trans;
    }
    P.get_floe_group().set_initial_window(new_window);
    P.make_input_file();

    return 0;
}