#include "../tests/catch.hpp"
#include <iostream>
#include "floe/geometry/geometries/point.hpp"
#include "floe/geometry/arithmetic/arithmetic.hpp"
#include "floe/geometry/arithmetic/dot_product.hpp"
#include "floe/geometry/arithmetic/point_operators.hpp"


TEST_CASE( "Test point arithmetic", "[geometry]" ) {

    using namespace std;
    using real = double;
    using point_type = floe::geometry::Point<real>;

    point_type P{1,1};
    // cout << norm2(P) << endl;
    REQUIRE(norm2(P) == sqrt(2));
    // cout << floe::geometry::dot_product(P, direct_orthogonal(P)) << endl;
    REQUIRE(floe::geometry::dot_product(P, direct_orthogonal(P)) == 0);

    point_type P2{2,3};
    point_type P3{3,4};
    cout << P + P2 << endl;
    // REQUIRE(P + P2 == P3);
}