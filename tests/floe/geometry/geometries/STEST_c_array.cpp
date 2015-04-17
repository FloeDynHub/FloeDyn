#include <iostream>
#include <type_traits>
#include <typeinfo>

#include <boost/concept_check.hpp>
#include <iostream>

#include <boost/geometry/util/bare_type.hpp>

#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/adapted/c_array.hpp"
#include "floe/geometry/geometries/point.hpp"
#include "floe/geometry/geometries/concepts/simple_static_polygon_concept.hpp"

BOOST_GEOMETRY_REGISTER_C_ARRAY_CS(cs::cartesian)

int main()
{
    using real = double;

    real triangle[3][2] = { {0,0}, {0,1}, {1,0} };

    using namespace std;
    namespace bg = boost::geometry;
    namespace fg = floe::geometry;

    using Point = fg::Point<real>;

    BOOST_CONCEPT_ASSERT( (fg::concept::SimpleStaticPolygon<decltype(triangle)>) );
    BOOST_CONCEPT_ASSERT( (fg::concept::ConstSimpleStaticPolygon<decltype(triangle)>) );

    cout << bg::get<2,0>(triangle) << endl;

    Point triangle2[3] = { {0, 0}, {0, 1}, {1, 0} };
    
    BOOST_CONCEPT_ASSERT( (fg::concept::SimpleStaticPolygon<decltype(triangle2)>) );
    BOOST_CONCEPT_ASSERT( (fg::concept::ConstSimpleStaticPolygon<decltype(triangle2)>) );
    
    
    cout << bg::get<2,0>(triangle2) << endl;

    cout << bg::dsv(triangle2[1]) << endl;
    cout << bg::dsv(triangle[1])  << endl;

    Point* triangle3[3] = { &triangle2[1], &triangle2[2], &triangle2[0] };
    BOOST_CONCEPT_ASSERT( (fg::concept::SimpleStaticPolygon<decltype(triangle3)>) );
    BOOST_CONCEPT_ASSERT( (fg::concept::ConstSimpleStaticPolygon<decltype(triangle3)>) );

    cout << bg::get<1,0>(triangle3) << endl;

    bg::set<1,0>(triangle3, 2.);
    cout << bg::dsv(triangle2[2]) << endl;

    return 0;
};

