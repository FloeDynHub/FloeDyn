#ifdef PBC

#include "../product/config/config_floes.hpp"
#include "floe/collision/matlab/periodic_detector.hpp"

namespace fcm = floe::collision::matlab;
using namespace types;

template class fcm::PeriodicMatlabDetector<floe_group_type, topology_type>;
using periodic_contact_type = typename fcm::PeriodicMatlabDetector<floe_group_type, topology_type>::contact_type;
using periodic_prox_data_type = typename fcm::PeriodicMatlabDetector<floe_group_type, topology_type>::proximity_data_type;
template class fcm::MatlabDetector<floe_group_type, periodic_prox_data_type, periodic_contact_type>;

#endif