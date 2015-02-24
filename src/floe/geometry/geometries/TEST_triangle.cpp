#include <iostream>

#include <boost/concept_check.hpp>

#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/adapted/c_array.hpp"
#include "floe/geometry/geometries/point.hpp"
#include "floe/geometry/geometries/triangle.hpp"
#include "floe/geometry/core/point_type.hpp"
#include "floe/geometry/core/coordinate_type.hpp"

int main() {
    using namespace std;
    namespace bg = boost::geometry;
    namespace fg = floe::geometry;

    using real = double;
    using Point = fg::Point<real>;

    Point triangle[3] = { {0,0}, {0,1}, {1,0} };

    cout << bg::static_num_points<decltype(triangle)>::value << endl;
    BOOST_CONCEPT_ASSERT( (bg::concept::Triangle<decltype(triangle)>) );

    using Triangle = fg::Triangle<Point>;
    Triangle triangle2 = { {0,0}, {0,1}, {1,0} };

    cout << bg::static_num_points<decltype(triangle2)>::value << endl;
    BOOST_CONCEPT_ASSERT( (bg::concept::Triangle<decltype(triangle2)>) );

    typedef typename fg::point_type<Triangle>::type pt_type;
    typedef typename fg::coordinate_type<Triangle>::type c_type;

    return 0;
}
