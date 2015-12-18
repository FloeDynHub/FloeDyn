#include "../product/config/config_dynamics.hpp"
#include "floe/ope/dynamics_manager.hpp"

template class floe::ope::DynamicsManager<external_forces_type, floe_group_type>;
template class floe::ope::DynamicsManager<generator_external_forces_type, floe_group_type>;