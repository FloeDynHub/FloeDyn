#include "../tests/catch.hpp"
#include <iostream>
#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/point.hpp" 
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/algorithms/intersects.hpp>
#include <cmath>

TEST_CASE( "Test boost geometry intersects algorithm", "[geometry]" ) {
    
    using namespace std;
    using point_type = floe::geometry::Point<double>;
    using Polygon = boost::geometry::model::polygon<point_type, false, false>;

    Polygon P0, P1, P2;
    uint N = 30;
    for (size_t i = 0; i != N; ++i)
    {
        double frac = (double)(2 * i * M_PI)/N;
        point_type t1 = {1, 0}, t2 = {2.1, 0};
        P0.outer().push_back({cos(frac), sin(frac)});
        P1.outer().push_back({t1.x + cos(frac), t1.y + sin(frac)});
        P2.outer().push_back({t2.x + cos(frac), t2.y + sin(frac)});
    }

    bool I1 = boost::geometry::intersects(P0, P1);
    bool I2 = boost::geometry::intersects(P0, P2);
    REQUIRE( I1 );
    REQUIRE( !I2 );

}
