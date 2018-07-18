#include <iostream>
#include "floe/geometry/geometry.hpp"
#include <boost/concept/assert.hpp>
#include "floe/geometry/geometries/point.hpp"
#include "floe/geometry/geometries/polygon.hpp"
#include "floe/geometry/geometries/box.hpp"
#include "floe/geometry/geometries/referring_box.hpp"

int main() {

    using namespace std;
    namespace bg = boost::geometry;
    namespace fg = floe::geometry;

    using real = double;
    using Point = floe::geometry::Point<real>;
    using Polygon = floe::geometry::Polygon<Point>;
    using Box = floe::geometry::Box<Point>;
    using ReferringBox = fg::ReferringBox<Point>;

    BOOST_CONCEPT_ASSERT( (boost::geometry::concepts::Polygon<Polygon>) );
    BOOST_CONCEPT_ASSERT( (boost::geometry::concepts::ConstPolygon<Polygon>) );

    auto polygon = Polygon();
    polygon.boundary().push_back({0.0,0.0});
    polygon.boundary().push_back({0.0,1.0});
    polygon.boundary().push_back({1.0,0.0});
    polygon.boundary().push_back({0.0,0.0});

    cout << bg::area(polygon) << endl;
    cout << bg::perimeter(polygon) << endl;
    cout << bg::dsv(bg::return_centroid<Point>(polygon)) << endl;
    
    cout << bg::distance(polygon, Point({1.,1.0})) << endl;

    auto box = bg::return_envelope<Box>(polygon);
    cout << bg::dsv(box) << endl;
   
    /*
    Box rbox = bg::return_envelope<ReferringBox>(polygon);
    cout << bg::dsv(rbox) << endl;
    */

    using Polygon2 = floe::geometry::Polygon<Point,true,false>;

    BOOST_CONCEPT_ASSERT( (boost::geometry::concepts::Polygon<Polygon2>) );

    auto polygon2 = Polygon2();
    polygon2.boundary().push_back({0.0,0.0});
    polygon2.boundary().push_back({0.0,1.0});
    polygon2.boundary().push_back({1.0,0.0});

    cout << bg::area(polygon2) << endl;
    cout << bg::perimeter(polygon2) << endl;
    cout << bg::dsv(bg::return_centroid<Point>(polygon2)) << endl;
    
    cout << bg::distance(polygon2, Point({1.,1.0})) << endl;

    auto box2 = bg::return_envelope<Box>(polygon2);
    cout << bg::dsv(box2) << endl;

    return 0;
}
