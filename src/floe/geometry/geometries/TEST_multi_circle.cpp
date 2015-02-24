#include <iostream>

#include <boost/math/constants/constants.hpp>

#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/point.hpp"
#include "floe/geometry/geometries/multi_circle.hpp"
#include "floe/geometry/geometries/circle.hpp"

#include "floe/geometry/frame/theta_frame.hpp"
#include "floe/geometry/frame/frame_transformers.hpp"

int main()
{
    namespace fg = floe::geometry;
    using real = double;
    using Point = fg::Point<real>;
    using Circle = fg::Circle<Point>;
    using MultiCircle = fg::MultiCircle<Circle>;

    const real pi = boost::math::constants::pi<real>();
    
    MultiCircle circles;
    circles.push_back( {{0.,0.}, 1.} );
    circles.push_back( {{1.,-1.}, 2.} );
    circles.push_back( {{2.,0.}, 1.} );

    using namespace std;
    cout << fg::dsv( circles ) << endl;

    using Frame = fg::frame::ThetaFrame<Point>;
    const Frame frame(1.0, 2.0, pi/4);

    MultiCircle circles2;
    fg::transform( circles, circles2, fg::frame::transformer(frame) );

    cout << fg::dsv( circles2 ) << endl;

    return 0;
}
