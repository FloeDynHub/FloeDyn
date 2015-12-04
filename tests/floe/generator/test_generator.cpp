#include "../tests/catch.hpp"
#include <iostream>
double DT_DEFAULT;
#include "../product/interrupt.hpp"
#include "../product/config_periodic.hpp"

#include "floe/generator/generator.hpp"



TEST_CASE( "Test floe generator", "[generator]" ) {

  problem_type P;
  floe::generator::Generator<problem_type> G(P);

  G.generate_floe_set(111, 0.7);

}
