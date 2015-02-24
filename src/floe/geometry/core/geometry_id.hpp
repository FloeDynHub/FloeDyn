/*!
 * \file floe/geometry/core/geometry_id.hpp
 * \brief Associate an unique id to each geometry. 
 * \author Roland Denis
 * \see boost/geometry/core/geometry_id.hpp
 */

#ifndef FLOE_GEOMETRY_CORE_GEOMETRY_ID_HPP
#define FLOE_GEOMETRY_CORE_GEOMETRY_ID_HPP

#include <boost/geometry/core/geometry_id.hpp>
#include "floe/geometry/core/tags.hpp"

namespace boost { namespace geometry
{

namespace core_dispatch
{

template <>
struct geometry_id<circle_tag> : boost::mpl::int_<192> {};

} // namespace core_dispatch

}} // namespace boost::geometry

#endif // FLOE_GEOMETRY_CORE_GEOMETRY_ID_HPP

