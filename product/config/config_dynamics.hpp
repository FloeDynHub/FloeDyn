#ifndef PRODUCT_CONFIG_CONFIG_DYNAMICS_HPP
#define PRODUCT_CONFIG_CONFIG_DYNAMICS_HPP

#include "../product/config/config_floes.hpp"

#include "floe/ope/physical_data.hpp"
#include "floe/ope/explicit_physical_data.hpp"
#include "floe/ope/external_forces.hpp"
#include "floe/ope/dynamics_manager.h"
#include "floe/domain/domain.hpp"

#ifdef PBC
#include "floe/ope/periodic_dynamics_manager.h"
#endif


using physical_data_type = floe::ope::PhysicalData<point_type>;
using external_forces_type = floe::ope::ExternalForces<floe_type, physical_data_type>;
using generator_physical_data_type = floe::ope::ExplicitPhysicalData<point_type>;
using generator_external_forces_type = floe::ope::ExternalForces<floe_type, generator_physical_data_type>;
using generator_dynamics_manager_type = floe::ope::DynamicsManager<generator_external_forces_type, floe_group_type>;

using domain_type = floe::domain::Domain<value_type>;

#ifdef PBC
using dynamics_manager_type = floe::ope::PeriodicDynamicsManager<external_forces_type, floe_group_type, topology_type>;
#else // Free boundary conditions types
using dynamics_manager_type = floe::ope::DynamicsManager<external_forces_type, floe_group_type>;
#endif


#endif // PRODUCT_CONFIG_CONFIG_DYNAMICS_HPP