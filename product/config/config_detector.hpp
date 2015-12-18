#ifndef PRODUCT_CONFIG_CONFIG_DETECTOR_HPP
#define PRODUCT_CONFIG_CONFIG_DETECTOR_HPP

#include "../product/config/config_floes.hpp"

#include "floe/dynamics/proximity_detector.hpp"
#include "floe/collision/matlab/detector.h"

#ifdef PBC
#include "floe/collision/matlab/periodic_detector.h"
#endif

using proximity_detector_type = floe::dynamics::ProximityDetector<
    floe::collision::matlab::MatlabDetector<floe_type>
>;
#ifdef PBC // Periodic boundary conditions types
using periodic_proximity_detector_type = floe::dynamics::ProximityDetector<
    floe::collision::matlab::PeriodicMatlabDetector<floe_type, topology_type>
>;
#endif


#endif // PRODUCT_CONFIG_CONFIG_DETECTOR_HPP