#include "../tests/catch.hpp"
#include <iostream>
#include "floe/geometry/geometries/point.hpp" 
#include "floe/dynamics/physical_data.hpp"


TEST_CASE( "Test physical data manager", "[ope]" ) {
    
    using namespace floe::dynamics;
    using namespace std;
    using point_type = floe::geometry::Point<double>;

    std::string file_name = "tests/floe/io/matlab/DataTopaz01.mat";

    double Time = 0;
    PhysicalData<point_type> P{Time};

    P.load_matlab_topaz_data(file_name);

    // for (std::size_t i = 0; i!= 10; ++i)
    //     { Time = 3600 * i; cout << P.water_speed() << " ";}
    // cout << endl;

    // for (std::size_t i = 0; i!= 60; ++i)
    //     { Time = 3600 * i; cout << P.water_speed() << " ";}
    int j = 40;
    Time = 3600 + 60 * j;
    auto pj = P.water_speed();
    Time = 3600;
    auto p1 = P.water_speed();
    Time = 3600 * 2;
    auto p2 = P.water_speed();
    REQUIRE(atan2(pj.y, pj.x) == ( 1 - (double)j/60) * atan2(p1.y, p1.x) + (double)j/60 * atan2(p2.y, p2.x));
    REQUIRE(norm2(pj) == ( 1 - (double)j/60) * norm2(p1) + (double)j/60 * norm2(p2));


    // Time = 165600.3098123;
    // for (int i = 0; i < 10; ++i)
    // {
    //     Time += 36.211098;
    //     cout << " W " << norm2(P.water_speed()) << " A " << norm2(P.air_speed()) << endl;
    // }

    double max_w = 0;
    int t0 = 1860, t1 = 1870;
    for (int i = 0; i < 5000; ++i)
        {Time = 60*i; max_w = std::max(max_w, norm2(P.water_speed())); if (norm2(P.water_speed()) > 5) cout << " ++ " << Time;}
    cout << " MAX W " << max_w;

    double sum_w = 0;
    for (int i = t0; i < t1; ++i)
        {Time = 60*i; sum_w += norm2(P.water_speed());}
    cout << " MOY W " << sum_w / (t1 - t0);

    double max_a = 0;
    for (int i = t0; i < t1; ++i)
        {Time = 60*i; max_a = std::max(max_a, norm2(P.air_speed()));}

    cout << " MAX A " << max_a;



}
