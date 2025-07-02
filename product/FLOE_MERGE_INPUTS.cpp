#include <iostream>
#include <cassert>
#include "../product/config/config.hpp"

/*
Create an hdf5 input file basically merging 2 input files
*/


int main( int argc, char* argv[] )
{   
    using namespace std;
    using namespace types;

    if ( argc < 3 )
    {
        cout << "Usage: " << argv[0] << " <input_file_name_1> <input_file_name_2>" << endl;
        return 1;
    }

    std::string input_list_floe_filename_1 = argv[1];
    std::string input_list_floe_filename_2 = argv[2];

    problem_type P;
    P.load_config(input_list_floe_filename_1);
    P.load_config(input_list_floe_filename_2);
    P.get_floe_group().set_initial_window(P.get_floe_group().bounding_window()); // TODO proper margin
    P.make_input_file();
    return 0;
}