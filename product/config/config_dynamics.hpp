#ifndef PRODUCT_CONFIG_CONFIG_DYNAMICS_HPP
#define PRODUCT_CONFIG_CONFIG_DYNAMICS_HPP

#include "../product/config/config_floes.hpp"

#include "floe/dynamics/physical_data.hpp"
#include "floe/dynamics/explicit_physical_data.hpp"
#include "floe/dynamics/external_forces.hpp"
#include "floe/dynamics/dynamics_manager.h"
#include "floe/domain/domain.hpp"

#ifdef PBC
#include "floe/dynamics/periodic_dynamics_manager.h"
#endif

namespace types {

// using physical_data_type = floe::dynamics::PhysicalData<point_type>;
using physical_data_type = floe::dynamics::ExplicitPhysicalData<point_type>;
using external_forces_type = floe::dynamics::ExternalForces<floe_type, physical_data_type>;
using generator_physical_data_type = floe::dynamics::ExplicitPhysicalData<point_type>;
using generator_external_forces_type = floe::dynamics::ExternalForces<floe_type, generator_physical_data_type>;
using generator_dynamics_manager_type = floe::dynamics::DynamicsManager<generator_external_forces_type, floe_group_type>;
// if physical_data_type == generator_physical_data_type
#define SAME_PHYSICAL_DATA_FOR_GENERATOR

using domain_type = floe::domain::Domain<value_type>;

#ifdef PBC
using dynamics_manager_type = floe::dynamics::PeriodicDynamicsManager<external_forces_type, floe_group_type, topology_type>;
#else // Free boundary conditions types
using dynamics_manager_type = floe::dynamics::DynamicsManager<external_forces_type, floe_group_type>;
#endif

} // namespace types

#endif // PRODUCT_CONFIG_CONFIG_DYNAMICS_HPP