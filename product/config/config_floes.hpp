#ifndef PRODUCT_CONFIG_CONFIG_FLOES_HPP
#define PRODUCT_CONFIG_CONFIG_FLOES_HPP

#include "../product/config/config_base.hpp"

#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"
#include "floe/variable/floe_group.hpp"

namespace ff = floe::floes;

using floe_type = ff::KinematicFloe<ff::StaticFloe<value_type, point_type>>;
using floe_group_type = floe::variable::FloeGroup<floe_type>;


#endif // PRODUCT_CONFIG_CONFIG_FLOES_HPP