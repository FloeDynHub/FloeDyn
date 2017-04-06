#include "../product/config/config_dynamics.hpp"
#include "floe/dynamics/dynamics_manager.hpp"
#include <type_traits>


template class floe::dynamics::DynamicsManager<types::external_forces_type, types::floe_group_type>;
#ifndef SAME_PHYSICAL_DATA_FOR_GENERATOR
template class floe::dynamics::DynamicsManager<types::generator_external_forces_type, types::floe_group_type>;
#endif