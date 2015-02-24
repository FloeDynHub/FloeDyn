/*!
 * \file floe/geometry/frame/TEST_theta_frame.cpp
 * \brief TEST of ThetaFrame class.
 * \author Roland Denis
 */

#include <iostream>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/point.hpp"
#include "floe/geometry/frame/theta_frame.hpp"
#include "floe/geometry/frame/frame_transformers.hpp"
#include "floe/geometry/geometries/circle.hpp"

int main()
{
    using namespace std;
    namespace fg = floe::geometry;

    using real = double;
    using Point = fg::Point<real>;
    using Frame = fg::frame::ThetaFrame<Point>;
    using Circle = fg::Circle<Point>;

    const real pi = boost::math::constants::pi<real>();

    const Frame frame(1.0, 2.0, pi/4);

    cout << "u = " << fg::dsv( frame.u() ) << " ; v = " << fg::dsv( frame.v() ) << endl;

    const auto strategy = fg::frame::transformer(frame);
    cout << "m = " << strategy.matrix() << endl;

    const Frame frame2(0., 1., pi/2);
    const auto strategy2 = fg::frame::transformer(frame,frame2);
    cout << "m2 = " << strategy2.matrix() << endl;

    Point pt;
    fg::transform( Point{1.,0.}, pt, strategy2 );
    cout << "pt = " << fg::dsv(pt) << endl;

    const Point pt2 = {1., 1.};
    const Circle circle{pt2, 2.};
    cout << fg::dsv(circle) << endl;
    Circle circle2;
    fg::transform( circle, circle2, strategy2 );
    cout << fg::dsv(circle2) << endl;

    return 0;
}
