#ifndef PRODUCT_CONFIG_CONFIG_DETECTOR_HPP
#define PRODUCT_CONFIG_CONFIG_DETECTOR_HPP

#include "../product/config/config_floes.hpp"

#include "floe/collision/matlab/detector.h"

#ifdef PBC
#include "floe/collision/matlab/periodic_detector.h"
#endif

#ifdef MPIRUN
#include "floe/collision/matlab/mpi_detector.hpp"
#endif

namespace types {

#ifdef MPIRUN
	using proximity_detector_type = floe::collision::matlab::MPIMatlabDetector<floe_group_type>;
#else
	using proximity_detector_type = floe::collision::matlab::MatlabDetector<floe_group_type>;
#endif

#ifdef PBC // Periodic boundary conditions types
	using periodic_proximity_detector_type = floe::collision::matlab::PeriodicMatlabDetector<floe_group_type, topology_type>;
#endif

}

#endif // PRODUCT_CONFIG_CONFIG_DETECTOR_HPP