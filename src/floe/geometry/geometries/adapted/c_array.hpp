// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2008-2012 Bruno Lalande, Paris, France.
// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2009-2012 Mateusz Loskot, London, UK.
// Copyright (c) 2014-     Roland Denis, Grenoble, France :p

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#ifndef FLOE_GEOMETRY_GEOMETRIES_ADAPTED_C_ARRAY_HPP
#define FLOE_GEOMETRY_GEOMETRIES_ADAPTED_C_ARRAY_HPP

// Check is original c_array.hpp has already been included
#ifdef BOOST_GEOMETRY_GEOMETRIES_ADAPTED_C_ARRAY_HPP
#error "boost/geometry/geometries/adapted/c_array.hpp must no be included."
#endif

#define BOOST_GEOMETRY_GEOMETRIES_ADAPTED_C_ARRAY_HPP // This file is incompatible with original from boost::geometry

#include <cstddef>

#include <boost/mpl/int.hpp>
#include <boost/type_traits/is_arithmetic.hpp>

#include "floe/geometry/core/access.hpp"
#include "floe/geometry/core/tags.hpp"
#include "floe/geometry/core/static_num_points.hpp"

#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/core/coordinate_dimension.hpp>
#include <boost/geometry/core/coordinate_type.hpp>

namespace boost { namespace geometry
{


#ifndef DOXYGEN_NO_TRAITS_SPECIALIZATIONS
namespace traits
{


#ifndef DOXYGEN_NO_DETAIL
namespace detail
{

// Detection of array of points
template <typename>
struct nested_c_array_tag
{
    typedef geometry_not_recognized_tag type;
};

template <>
struct nested_c_array_tag<point_tag>
{
    typedef simple_static_polygon_tag type;
};


// Create class and specialization to indicate the tag
// for normal cases and the case that the type of the c-array is arithmetic
template <bool, typename>
struct c_array_tag
{
    typedef geometry_not_recognized_tag type;
};


template <typename T>
struct c_array_tag<true, T>
{
    typedef point_tag type;
};

template <typename T>
struct c_array_tag<false, T>
    : nested_c_array_tag<typename geometry::tag<T>::type>
{};

} // namespace detail
#endif // DOXYGEN_NO_DETAIL


// Assign the point-tag, preventing arrays of points getting a point-tag
template <typename CoordinateType, std::size_t DimensionCount>
struct tag<CoordinateType[DimensionCount]>
    : detail::c_array_tag<boost::is_arithmetic<CoordinateType>::value, CoordinateType> {};

template <typename CoordinateType, std::size_t DimensionCount>
struct coordinate_type<CoordinateType[DimensionCount]>
{
    typedef CoordinateType type;
};

template <typename CoordinateType, std::size_t DimensionCount>
struct dimension<CoordinateType[DimensionCount]>: boost::mpl::int_<DimensionCount> {};


template <typename CoordinateType, std::size_t DimensionCount, std::size_t Dimension>
struct access<CoordinateType[DimensionCount], Dimension>
{
    static inline CoordinateType get(CoordinateType const p[DimensionCount])
    {
        return p[Dimension];
    }

    static inline void set(CoordinateType p[DimensionCount],
        CoordinateType const& value)
    {
        p[Dimension] = value;
    }
};

template <typename PointType, std::size_t NumPoints, std::size_t Index, std::size_t Dimension>
struct indexed_access<PointType[NumPoints], Index, Dimension>
{
    typedef typename geometry::coordinate_type<PointType>::type coordinate_type;

    static inline 
    coordinate_type get( PointType const p[NumPoints])
    {
        return geometry::get<Dimension>( p[Index] );
    }

    static inline
    void set( PointType p[NumPoints], coordinate_type  const& value )
    {
        geometry::set<Dimension>( p[Index], value );
    }
};

template <typename PointType, std::size_t NumPoints>
struct point_type<PointType[NumPoints]>
{
    typedef PointType type;
};


template <
    typename CoordinateType, 
    std::size_t NumPoints
>
struct static_num_points<CoordinateType[NumPoints]>
    : boost::mpl::int_<NumPoints> 
{};

} // namespace traits
#endif // DOXYGEN_NO_TRAITS_SPECIALIZATIONS


}} // namespace boost::geometry


#define BOOST_GEOMETRY_REGISTER_C_ARRAY_CS(CoordinateSystem) \
    namespace boost { namespace geometry { namespace traits { \
    template <typename T, std::size_t N> \
    struct coordinate_system<T[N]> \
    { \
        typedef CoordinateSystem type; \
    }; \
    }}}


#endif // BOOST_GEOMETRY_GEOMETRIES_ADAPTED_C_ARRAY_HPP
