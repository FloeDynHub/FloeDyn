#ifndef PRODUCT_CONFIG_CONFIG_BASE_HPP
#define PRODUCT_CONFIG_CONFIG_BASE_HPP
#include "floe/geometry/geometries/point.hpp"
#ifdef PBC
#include "floe/topology/toric_topology.hpp"
#endif
namespace types {

using value_type = double;
using point_type = floe::geometry::Point<value_type>;

#ifdef PBC
using topology_type = floe::topology::ToricTopology<point_type>;
#endif

} // namespace types

#endif // PRODUCT_CONFIG_CONFIG_BASE_HPP