#include "../product/config/config_dynamics.hpp"
#include "../product/config/config_floes.hpp"
#include "floe/io/hdf5_manager.hpp"
#include <type_traits>

template class floe::io::HDF5Manager<types::floe_group_type, types::dynamics_manager_type>;
#if !defined( SAME_PHYSICAL_DATA_FOR_GENERATOR ) || defined( PBC )
template class floe::io::HDF5Manager<types::floe_group_type, types::generator_dynamics_manager_type>;
#endif

