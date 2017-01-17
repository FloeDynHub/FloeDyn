#include "../tests/catch.hpp"

#include <iostream>
#include "floe/topology/toric_topology.hpp"
#include "floe/geometry/geometry.hpp"

// #include <chrono>


TEST_CASE( "Test topology", "[topology]" ) {

    using namespace std;
    using namespace floe::topology;
    using real = double;
    using point_type = floe::geometry::Point<real>;

    real x0 = 0, x1 = 10, y0 = 0, y1 = 20; 
    ToricTopology<real, point_type> T{x0, x1, y0, y1};

    point_type pt{x1 + 1, y0 - 1};
    T.replace(pt);
    REQUIRE( equal_points( pt, point_type{1, y1 -1}) );

    auto g_list = T.ghosts( pt );
    cout << "ghosts for " << pt << " : ";
    REQUIRE( g_list.size() == 4);
    for (auto& pt : g_list)
        cout << pt;

    real dx = x1 - x0, dy = y1 - y0;
    REQUIRE( equal_points(g_list[0] - pt, {dx, 0}));
    REQUIRE( equal_points(g_list[1] - pt, {dx, dy}));
    REQUIRE( equal_points(g_list[2] - pt, {0, dy}));
    REQUIRE( equal_points(g_list[3] - pt, {-dx, dy}));


}