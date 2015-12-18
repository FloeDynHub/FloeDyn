#ifdef PBC

#include "../product/config/config_dynamics.hpp"
#include "floe/dynamics/periodic_dynamics_manager.hpp"

template class floe::dynamics::PeriodicDynamicsManager<external_forces_type, floe_group_type, topology_type>;

#endif