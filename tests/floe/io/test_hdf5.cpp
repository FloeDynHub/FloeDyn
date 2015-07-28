#include "../tests/catch.hpp"

#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"
#include "floe/variable/floe_group.hpp"

#include <iostream>
#include <string>

namespace ff = floe::floes;

TEST_CASE( "Test hdf5 output", "[io]" )
{

    using namespace floe::variable;
    using namespace std;
    using floe_type = ff::KinematicFloe<ff::StaticFloe<double>>;
    using point_type = typename floe_type::point_type;
    using value_type = typename floe_type::value_type;
    using mesh_type = floe_type::mesh_type;
    FloeGroup<floe_type> F;
    // Import floes from Matlab configuration
    std::string mat_file_name = "tests/floe/io/matlab/r1day_set_up_250sm_sz_60_list_so_350_str.mat";
    // std::string mat_file_name = "tests/floe/io/matlab/config_q2_str.mat";
    F.load_matlab_config(mat_file_name);

    F.recover_states_from_file("io/out.h5", 5.23);

    // for (int i=0; i < 100; ++i)
    //     F.out_hdf5(i);
}


