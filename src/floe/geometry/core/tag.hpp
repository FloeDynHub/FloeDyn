/*!
 * \file floe/geometry/core/tag.hpp
 * \brief Utilities to check tag of geometries. Simply import boost version.
 * \author Roland Denis
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/geometry/reference/core/tag.html">boost/geometry/core/tag.hpp</a>
 */

#ifndef FLOE_GEOMETRY_CORE_TAG_HPP
#define FLOE_GEOMETRY_CORE_TAG_HPP

#include "floe/geometry/core/tags.hpp"

#include <boost/geometry/core/tag.hpp>

namespace floe { namespace geometry {

using boost::geometry::tag;

}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_CORE_TAG_HPP
