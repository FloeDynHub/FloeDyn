#include <iostream>

#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/point.hpp"
#include "floe/geometry/geometries/circle.hpp"
#include "floe/geometry/geometries/box.hpp"
#include "floe/geometry/geometries/segment.hpp"

int main() {

    using namespace std;
    namespace bg = boost::geometry;
    namespace fg = floe::geometry;

    using real = double;
    using Point = fg::Point<real>;
    using Circle = fg::Circle<Point>;
    using Box = fg::Box<Point>;

    BOOST_CONCEPT_ASSERT( (boost::geometry::concept::Circle<Circle>) );
    
    Circle circle{ {1.0,0.0}, 2. };

    cout << fg::area(&circle) << endl;
    cout << bg::dsv( fg::return_envelope<Box>(circle) ) << endl;
    cout << bg::dsv( fg::return_centroid<Point>(circle) ) << endl;
    cout << get_radius( circle ) << endl;
    cout << get<0>(&circle) << endl;
    
    cout << bg::dsv( fg::return_envelope<Box>( fg::return_buffer<Circle>( circle, 1 ) ) ) << endl;
    cout << fg::dsv( fg::return_buffer<Circle>( circle, 1 ) ) << endl;

    return 0;
}
