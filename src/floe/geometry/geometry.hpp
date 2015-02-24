/*!
 * \file floe/geometry/geometry.hpp
 * \brief To includes main files for manipulating geometries
 * \see boost/geometry/geometry.hpp
 */

#ifndef FLOE_GEOMETRY_HPP
#define FLOE_GEOMETRY_HPP

#include "floe/geometry/core/tag.hpp"
#include "floe/geometry/core/tags.hpp"

// Core algorithms
#include "floe/geometry/core/radius.hpp"
#include "floe/geometry/core/access.hpp"
#include "floe/geometry/core/static_num_points.hpp"
#include "floe/geometry/core/exterior_ring.hpp"
#include "floe/geometry/core/ring_type.hpp"
#include "floe/geometry/core/coordinate_dimension.hpp"
#include "floe/geometry/core/coordinate_type.hpp"
#include "floe/geometry/core/cells_type.hpp"
#include "floe/geometry/core/cells.hpp"
#include "floe/geometry/core/point_type.hpp"
#include "floe/geometry/core/point_order.hpp"

#include "floe/geometry/arithmetic/arithmetic.hpp"
#include "floe/geometry/arithmetic/dot_product.hpp"

// Algorithms
#include "floe/geometry/algorithms/area.hpp"
#include "floe/geometry/algorithms/envelope.hpp"
#include "floe/geometry/algorithms/centroid.hpp"
#include "floe/geometry/algorithms/buffer.hpp"
#include "floe/geometry/algorithms/transform.hpp"
#include "floe/geometry/algorithms/assign.hpp"
#include "floe/geometry/algorithms/num_points.hpp"
#include "floe/geometry/algorithms/distance.hpp"
#include "floe/geometry/algorithms/comparable_distance.hpp"
#include "floe/geometry/algorithms/disjoint.hpp"


// Concepts (check includes all concepts)
#include "floe/geometry/geometries/concepts/check.hpp"

// Views
#include "floe/geometry/views/center_view.hpp"

// IO
#include "floe/geometry/io/dsv/write.hpp"

// boost::geometry
#include <boost/geometry.hpp>

#endif // FLOE_GEOMETRY_HPP

/*!
 * \namespace floe::geometry
 * \brief Contains all stuffs about geometries.
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/index.html">boost/geometry.hpp</a>
 *
 * \namespace boost::geometry
 * \brief Contains all stuffs about geometries.
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/index.html">boost/geometry.hpp</a>
 *
 * \namespace floe::geometry::detail
 * \brief Implementation details about geometries tools.
 *
 * \namespace boost::geometry::detail
 * \brief Implementation details about geometries tools.
 *
 * \namespace boost::geometry::dispatch
 * \brief Dispatch functions for geometry utilities.
 *
 * \namespace boost::geometry::traits
 * \brief Type traits about geometries.
 */
