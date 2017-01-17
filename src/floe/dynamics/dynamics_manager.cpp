#include "../product/config/config_dynamics.hpp"
#include "floe/dynamics/dynamics_manager.hpp"
#include <type_traits>


template class floe::dynamics::DynamicsManager<types::external_forces_type, types::floe_group_type>;
// constexpr if (!std::is_same<types::external_forces_type, types::generator_external_forces_type>::value)
template class floe::dynamics::DynamicsManager<types::generator_external_forces_type, types::floe_group_type>;

// Trying to conditionally explicit instantiate a type :
// using foo = typename std::conditional<
// 	std::is_same<types::external_forces_type, types::generator_external_forces_type>::value,
// 	void,
// 	floe::dynamics::DynamicsManager<types::generator_external_forces_type, types::floe_group_type>
// >::type;
