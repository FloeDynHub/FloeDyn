#include "../tests/catch.hpp"
#include <iostream>
#include "floe/geometry/geometries/point.hpp" 
#include "floe/ope/physical_data.hpp"


TEST_CASE( "Test physical data manager", "[ope]" ) {
    
    using namespace floe::ope;
    using namespace std;
    using point_type = floe::geometry::Point<double>;

    std::string file_name = "tests/floe/io/matlab/DataTopaz01.mat";

    PhysicalData<point_type> P;

    P.import_topaz_mat_file(file_name);

    for (std::size_t i = 0; i!= 10; ++i)
        cout << P.water_speed(3600 * i, {0,0}) << " ";
    cout << endl;

    for (std::size_t i = 0; i!= 60; ++i)
        cout << P.water_speed(60 * i, {0,0}) << " ";

    int j = 40;
    auto pj = P.water_speed(3600 + 60 * j);
    auto p1 = P.water_speed(3600);
    auto p2 = P.water_speed(3600 * 2);
    REQUIRE(atan2(pj.y, pj.x) == ( 1 - (double)j/60) * atan2(p1.y, p1.x) + (double)j/60 * atan2(p2.y, p2.x));
    REQUIRE(norm2(pj) == ( 1 - (double)j/60) * norm2(p1) + (double)j/60 * norm2(p2));

}
