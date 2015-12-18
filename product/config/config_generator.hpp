#ifndef PRODUCT_CONFIG_CONFIG_GENERATOR_HPP
#define PRODUCT_CONFIG_CONFIG_GENERATOR_HPP

#include "../product/config/config_floes.hpp"
#include "../product/config/config_detector.hpp"
#include "../product/config/config_collision.hpp"
#include "../product/config/config_dynamics.hpp"
#include "floe/problem/problem.hpp"
#include "floe/generator/generator.h"

namespace types {

using generator_problem_type = floe::problem::Problem<
    floe_group_type,
    proximity_detector_type,
    generator_collision_manager_type,
    generator_dynamics_manager_type,
    domain_type
>;
using generator_type = floe::generator::Generator<generator_problem_type>;

}

#endif // PRODUCT_CONFIG_CONFIG_GENERATOR_HPP