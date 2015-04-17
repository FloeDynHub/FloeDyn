/*!
 * \file floe/state/TEST_space_time_state.cpp
 * \brief Test for SpaceTimeState class.
 * \author Roland Denis
 */

#include <iostream>
#include "floe/geometry/arithmetic/point_operators.hpp"
#include "floe/geometry/geometries/point.hpp"
#include "space_time_state.hpp"

int main()
{
    using namespace std;
    using real = double;
    using Point = floe::geometry::Point<real>;
    using State = floe::state::SpaceTimeState<Point,real>;

    State state;
    state.pos.x = 0.; state.pos.y = 1.;
    state.theta = 1.;
    state.speed.x = 2.; state.speed.y = 2.;
    state.rot = 90.;
    cout << state << endl;
    cout << 2*state << endl;
    cout << (state + 2*state) << endl;
    return 0;
}
