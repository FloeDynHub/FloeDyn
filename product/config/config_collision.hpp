#ifndef PRODUCT_CONFIG_CONFIG_COLLISION_HPP
#define PRODUCT_CONFIG_CONFIG_COLLISION_HPP

#include "../product/config/config_base.hpp"


#include "floe/lcp/LCP_manager.hpp"
#include "floe/lcp/solver/LCP_solver.hpp"
#include "floe/lcp/solver/generator_LCP_solver.hpp"
#include "floe/collision/collision_manager.hpp"

using solver_type = floe::lcp::solver::LCPSolver<value_type>;
using manager_h_type = floe::lcp::LCPManager<solver_type>;
using collision_manager_type = floe::collision::CollisionManager<manager_h_type>;
using generator_solver_type = floe::lcp::solver::GeneratorLCPSolver<value_type>;
using generator_manager_h_type = floe::lcp::LCPManager<generator_solver_type>;
using generator_collision_manager_type = floe::collision::CollisionManager<generator_manager_h_type>;


#endif // PRODUCT_CONFIG_CONFIG_COLLISION_HPP