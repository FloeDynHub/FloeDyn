#include "../tests/catch.hpp"

#include "../product/config/config_dynamics.hpp"
#include "floe/io/hdf5_manager.h"

#include <iostream>
#include <string>

namespace ff = floe::floes;

TEST_CASE( "Test hdf5 io", "[io]" )
{
    floe_group_type F;
    double t = 0;
    dynamics_manager_type D(t);
    // Import floes from Matlab configuration
    std::string mat_file_name = "tests/floe/io/matlab/r1day_set_up_250sm_sz_60_list_so_350_str.mat";
    // std::string mat_file_name = "tests/floe/io/matlab/config_q2_str.mat";
    F.load_matlab_config(mat_file_name);

    auto hdf_mgr = floe::io::HDF5Manager<floe_group_type, dynamics_manager_type>(F);

    // F.recover_states_from_file("io/out.h5", 5.23);

    for (int i=0; i < 100; ++i)
        hdf_mgr.save_step(i, D);
}


