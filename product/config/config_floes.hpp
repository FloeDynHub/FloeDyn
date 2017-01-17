#ifndef PRODUCT_CONFIG_CONFIG_FLOES_HPP
#define PRODUCT_CONFIG_CONFIG_FLOES_HPP

#include "../product/config/config_base.hpp"

#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"
#include "floe/floes/floe_group.hpp"

#ifdef MPIRUN
#include "floe/floes/partial_floe_group.hpp"
#include "floe/floes/identifiable_mixin.hpp"
#endif
// #include "floe/floes/partial_floe_group.hpp"

namespace ff = floe::floes;

namespace types { // testings


#ifdef MPIRUN
using floe_type = ff::Identifiable<std::size_t, ff::KinematicFloe<ff::StaticFloe<value_type, point_type>>>;
using floe_group_type = floe::floes::PartialFloeGroup<floe_type>;
#else
using floe_type = ff::KinematicFloe<ff::StaticFloe<value_type, point_type>>;
using floe_group_type = floe::floes::FloeGroup<floe_type>;
// using floe_group_type = floe::floes::PartialFloeGroup<floe_type>;
#endif


} // namespace types

#endif // PRODUCT_CONFIG_CONFIG_FLOES_HPP