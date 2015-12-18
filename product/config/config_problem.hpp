#ifndef PRODUCT_CONFIG_CONFIG_PROBLEM_HPP
#define PRODUCT_CONFIG_CONFIG_PROBLEM_HPP

#include "../product/config/config_floes.hpp"
#include "../product/config/config_detector.hpp"
#include "../product/config/config_collision.hpp"
#include "../product/config/config_dynamics.hpp"
#include "floe/problem/problem.hpp"

#ifdef PBC
#include "floe/problem/periodic_problem.hpp"
#endif

namespace types {

using namespace floe::problem;

#ifdef PBC // Periodic boundary conditions types
using problem_type = PeriodicProblem<
    floe_group_type,
    periodic_proximity_detector_type,
    collision_manager_type,
    dynamics_manager_type,
    domain_type,
    topology_type
>;
#else // Free boundary conditions types
using problem_type = Problem<
    floe_group_type,
    proximity_detector_type,
    collision_manager_type,
    dynamics_manager_type,
    domain_type
>;
#endif

}

#endif // PRODUCT_CONFIG_CONFIG_PROBLEM_HPP