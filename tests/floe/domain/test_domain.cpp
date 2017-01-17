#include "../tests/catch.hpp"
#include <iostream>
#include "floe/domain/domain.hpp"


TEST_CASE( "Test Domain", "[doman]" ) {
    using namespace floe::domain;
    Domain D;
    REQUIRE( D.time() == 0 );
    REQUIRE( D.time_step() == 1 );
    D.set_time_step(0.1);
    REQUIRE( D.time_step() == 0.1 );
}