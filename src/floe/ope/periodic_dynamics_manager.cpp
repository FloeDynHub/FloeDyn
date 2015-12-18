#ifdef PBC

#include "../product/config/config_dynamics.hpp"
#include "floe/ope/periodic_dynamics_manager.hpp"

template class floe::ope::PeriodicDynamicsManager<external_forces_type, floe_group_type, topology_type>;

#endif