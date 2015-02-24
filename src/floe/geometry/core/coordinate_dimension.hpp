/*!
 * \file floe/geometry/core/coordinate_dimension.hpp
 * \brief Dimension getter and comparators from boost.
 * \author Roland Denis
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/geometry/reference/core/dimension.html">boost/geometry/core/dimension.hpp</a>
 */

#ifndef FLOE_GEOMETRY_CORE_COORDINATE_DIMENSION_HPP
#define FLOE_GEOMETRY_CORE_COORDINATE_DIMENSION_HPP

#include <boost/geometry/core/coordinate_dimension.hpp>

namespace floe { namespace geometry
{
    using boost::geometry::dimension;
    using boost::geometry::assert_dimension;
    using boost::geometry::assert_dimension_less_equal;
    using boost::geometry::assert_dimension_greater_equal;
    using boost::geometry::assert_dimension_equal;

}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_CORE_COORDINATE_DIMENSION_HPP

