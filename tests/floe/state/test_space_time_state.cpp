#include "../tests/catch.hpp"

#include <iostream>
#include "floe/geometry/geometries/point.hpp"
#include "floe/state/space_time_state.hpp"

TEST_CASE( "Test space time state", "[state]" )
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

    State state2{{1, 2}, 3, {4, 5}, 6};
    cout << state2 << endl;
    REQUIRE(state2.pos.x == 1);
    REQUIRE(state2.pos.y == 2);
    REQUIRE(state2.theta == 3);
    REQUIRE(state2.speed.x == 4);
    REQUIRE(state2.speed.y == 5);
    REQUIRE(state2.rot == 6);
}
