#include "../tests/catch.hpp"

#include <iostream>

#include "floe/floes/ghost_floe.hpp"
#include "floe/floes/floe_group.hpp"

// #include <chrono>
#include <algorithm>


TEST_CASE( "Test ghost floe", "[floes]" ) {

    namespace ff = floe::floes;
    using namespace std;

    using floe_type = ff::KinematicFloe<ff::StaticFloe<double>>;
    using point_type = typename floe_type::point_type;

    floe::floes::FloeGroup<floe_type> F;

    // Import floes from Matlab configuration
    std::string mat_file_name = "tests/floe/io/matlab/1floe_rot1.mat";
    F.load_matlab_config(mat_file_name);

    floe_type& floe = F.get_floes()[0];
    point_type translation = point_type{3.5, 7};
    ff::GhostFloe<floe_type> ghost( floe, translation);

    REQUIRE(equal_points(ghost.frame().center(), floe.frame().center() + translation));
    REQUIRE(equal_points(ghost.state().pos, floe.state().pos + translation));
    
    REQUIRE(ghost.geometry().outer().size() == floe.geometry().outer().size());

    vector<bool> eq(ghost.geometry().outer().size());
    for (std::size_t i = 0; i != eq.size(); ++i)
        eq[i] = equal_points(ghost.geometry().outer()[i], floe.geometry().outer()[i] + translation);

    REQUIRE( std::all_of(eq.begin(), eq.end(), [](bool b){ return b; }));

    // auto t_start = chrono::high_resolution_clock::now();
    // for (int i = 0; i < 1000 ; ++i)
    // {
    //     cout << ghost.geometry().outer()[0].x;
    // }
    // auto t_end = chrono::high_resolution_clock::now();
    // cout << "get ghost geometry : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;

}