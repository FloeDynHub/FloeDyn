#include "../product/config/config_floes.hpp"
#include "../product/config/config_dynamics.hpp"
#include "floe/io/hdf5_manager.hpp"

template class floe::io::HDF5Manager<floe_group_type, dynamics_manager_type>;
template class floe::io::HDF5Manager<floe_group_type, generator_dynamics_manager_type>;