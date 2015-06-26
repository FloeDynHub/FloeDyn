#ifndef FLOE_GEOMETRY_ARITHMETIC_ARITHMETIC_HPP
#define FLOE_GEOMETRY_ARITHMETIC_ARITHMETIC_HPP

#include <boost/geometry/arithmetic/arithmetic.hpp>
#include <math.h>  

namespace floe { namespace geometry
{

using boost::geometry::add_value;
using boost::geometry::add_point;
using boost::geometry::subtract_value;
using boost::geometry::subtract_point;
using boost::geometry::multiply_value;
using boost::geometry::multiply_point;
using boost::geometry::divide_value;
using boost::geometry::divide_point;
using boost::geometry::assign_value;
using boost::geometry::assign_point;


template<typename TPoint>
TPoint direct_orthogonal(TPoint const& point)
{
    return TPoint{-point.y, point.x};
}

template<typename TPoint>
auto cross_product_value(TPoint const& point1, TPoint const& point2)
-> decltype(point1.x)
{
    return point1.x * point2.y -point1.y * point2.x;
}



}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_ARITHMETIC_ARITHMETIC_HPP

