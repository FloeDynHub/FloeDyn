// test file for Xcode debugging
#include <iostream>
#include "../tests/floe/config.hpp"


namespace ff = floe::floes;

int main( int argc, char* argv[] )
{
    std::string mat_file_name = "/Users/Serge/Projects/FloeCpp/tests/floe/io/matlab/r1day_set_up_250sm_sz_60_list_so_350_str.mat";
    problem_type P;
    P.create_h();
    P.load_matlab_config(mat_file_name);
    P.solve();

    return 0;
}
