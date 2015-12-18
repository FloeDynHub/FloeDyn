#include "../product/config/config_dynamics.hpp"
#include "floe/dynamics/dynamics_manager.hpp"

template class floe::dynamics::DynamicsManager<external_forces_type, floe_group_type>;
template class floe::dynamics::DynamicsManager<generator_external_forces_type, floe_group_type>;