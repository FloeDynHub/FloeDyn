#include "../tests/catch.hpp"
#include <iostream>

#include "floe/variable/floe_group.hpp"
#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"

namespace ff = floe::floes;

TEST_CASE( "Test Floes", "[variable]" ) {

    using namespace floe::variable;
    using namespace std;

    using floe_type = ff::KinematicFloe<ff::StaticFloe<double>>;
    using mesh_type = floe_type::mesh_type;

    // using floe_alg_type = Floe_alg;
    // using floe_h_type = Floe_h<mesh_type, floe_alg_type>;
    // using floe_group_h_type = FloeGroup_h<floe_h_type>;

    // Create empty Floe Group
    // FloeGroup<floe_type, floe_group_h_type> F;
    FloeGroup<floe_type> F;
    REQUIRE( F.m_list_floe.size() == 0 );

    // Import floes from Matlab configuration
    std::string mat_file_name = "tests/floe/io/matlab/r1day_set_up_250sm_sz_60_list_so_350_str.mat";
    F.load_matlab_config(mat_file_name);
    REQUIRE(F.m_list_floe.size() == 350);
    // is discrete group correctly instanciated ?
    REQUIRE(F.get_floe_group_h().m_list_floe_h.size() == 350);

}