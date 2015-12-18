#ifndef PRODUCT_CONFIG_CONFIG_COLLISION_HPP
#define PRODUCT_CONFIG_CONFIG_COLLISION_HPP

#include "../product/config/config_base.hpp"


#include "floe/ope/LCP_manager.hpp"
#include "floe/ope/LCP_solver.hpp"
#include "floe/ope/generator_LCP_solver.hpp"
#include "floe/ope/collision_manager.hpp"

using solver_type = floe::ope::LCPSolver<value_type>;
using manager_h_type = floe::ope::LCPManager<solver_type>;
using collision_manager_type = floe::ope::CollisionManager<manager_h_type>;
using generator_solver_type = floe::ope::GeneratorLCPSolver<value_type>;
using generator_manager_h_type = floe::ope::LCPManager<generator_solver_type>;
using generator_collision_manager_type = floe::ope::CollisionManager<generator_manager_h_type>;


#endif // PRODUCT_CONFIG_CONFIG_COLLISION_HPP