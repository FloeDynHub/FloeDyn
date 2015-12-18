#ifdef PBC

#include "../product/config/config_floes.hpp"
#include "floe/collision/matlab/periodic_detector.hpp"

namespace fcm = floe::collision::matlab;

template class fcm::PeriodicMatlabDetector<floe_type, topology_type>;
using periodic_contact_type = typename fcm::PeriodicMatlabDetector<floe_type, topology_type>::contact_type;
template class fcm::MatlabDetector<floe_type, periodic_contact_type>;

#endif