#include <iostream>
#include <boost/geometry.hpp>
#include <boost/concept/assert.hpp>
#include "floe/geometry/geometries/point.hpp"
#include "floe/geometry/arithmetic/point_operators.hpp"

int main () {
    using namespace std;
    using real = double;
    using Point = floe::geometry::Point<real>;

    BOOST_CONCEPT_ASSERT( (boost::geometry::concept::Point<Point>) );

    auto pt = Point(1.0, 0.0);
    get<1>(pt) = 2.0;

    cout << pt.x << " ; " << pt.y << endl;

    if (pt.y == 2.0) {
        cout << "OK" << endl;
    } else {
        cout << "KO" << endl;
    }

    const auto pt2 = Point(0.0, 1.0);
    cout << boost::geometry::distance(pt, pt2) << endl;

    cout << norm2(pt2 - pt) << endl;

    return 0;
}
