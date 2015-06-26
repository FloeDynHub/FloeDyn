#include "../tests/catch.hpp"

#include <iostream>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/point.hpp"
#include "floe/geometry/frame/theta_frame.hpp"
#include "floe/geometry/frame/frame_transformers.hpp"
#include "floe/geometry/geometries/circle.hpp"

#include <chrono>

TEST_CASE( "Test frame transformer", "[frame]" ) {

    using namespace floe::geometry;
    using namespace std;
    namespace trans = boost::geometry::strategy::transform;

    using point_type = Point<double>;
    using frame_type = frame::ThetaFrame<point_type>;
    using polygon = boost::geometry::model::polygon<point_type, false, false>;


    const auto transf = frame::transformer( frame_type{point_type{2,4}, 0.7} );

    polygon P{};
    polygon P1{};

    std::vector<point_type> V{
        point_type{0,0},
        point_type{1,0},
        point_type{1,1},
        point_type{0,1}
    };

    for (auto& point: V)
        P.outer().push_back(point);

    cout << dsv(P) << endl;

    double D1{distance(P.outer()[0], P.outer()[2])};
    transform( P, P1, transf );

    cout << dsv(P1) << endl;

    double D2{distance(P.outer()[0], P.outer()[2])};

    REQUIRE(D1 == D2);


    int N = 2000;
    polygon B{}, B2{};
    for (int i = 0; i < N ; ++i)
        B.outer().push_back(point_type{cos((double) i / N), sin((double) i / N)});


    // CHRONO TRANSFORM
    auto t_start = chrono::high_resolution_clock::now();
    transform( B, B2, frame::transformer( frame_type{point_type{1.23, 3.21}, 0} ) );
    auto t_end = chrono::high_resolution_clock::now();
    cout << "boost transform + copy : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;


    t_start = chrono::high_resolution_clock::now();
    transform( B, B, frame::transformer( frame_type{point_type{1.23, 3.21}, 3} ) );
    t_end = chrono::high_resolution_clock::now();
    cout << "transform : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;

    t_start = chrono::high_resolution_clock::now();
    B2 = B;
    t_end = chrono::high_resolution_clock::now();
    cout << "copy : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;

    t_start = chrono::high_resolution_clock::now();
    point_type T{1.23, 3.21};
    for (point_type& pt : B.outer())
        add_point(pt, T);
    t_end = chrono::high_resolution_clock::now();
    cout << "manual translate : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;

    t_start = chrono::high_resolution_clock::now();
    trans::translate_transformer<double, 2,2> translate(1.5, 1.5);
    boost::geometry::transform(B, B2, translate);
    t_end = chrono::high_resolution_clock::now();
    cout << "boost translate : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;
    
    t_start = chrono::high_resolution_clock::now();
    trans::rotate_transformer<boost::geometry::degree, double, 2, 2> rotate(90.0);
    boost::geometry::transform(B, B2, rotate);
    t_end = chrono::high_resolution_clock::now();
    cout << "boost rotate : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;
    // CHRONO TRANSFORM
}