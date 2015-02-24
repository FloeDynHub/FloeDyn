#ifndef FLOE_GEOMETRY_GEOMETRIES_TRIANGLE_HPP
#define FLOE_GEOMETRY_GEOMETRIES_TRIANGLE_HPP

#include <cstddef>
#include <array>

#include <boost/mpl/int.hpp>

#include "floe/geometry/core/tag.hpp"
#include "floe/geometry/core/access.hpp"
#include "floe/geometry/core/static_num_points.hpp"
#include "floe/geometry/core/coordinate_type.hpp"

namespace floe { namespace geometry
{

/*! Triangle model using std::array
 *
 * A triangle is a model of SimpleStaticPolygon composed of 3 points
 *
 * \tparam TPoint Point type
 * \tparam ClockWise true if the points are clockwise oriented
 */
template <
    typename TPoint,
    bool ClockWise = true
>
struct Triangle
    : std::array<TPoint, 3>
{
    Triangle( TPoint const& p1, TPoint const& p2, TPoint const& p3 )
        : std::template array<TPoint,3>{ { p1, p2, p3  }}
    {}
};

}} // namespace floe::geometry

// Type traits

// For easy use, locally
namespace {
    using floe::geometry::Triangle;
}

// Boost geometyr traits
namespace boost { namespace geometry { namespace traits
{

//!Tag
template <
    typename TPoint,
    bool ClockWise
>
struct tag< Triangle<TPoint,ClockWise> >
{
    typedef triangle_tag type;
};

template <
    typename TPoint,
    bool ClockWise
>
struct static_num_points< Triangle<TPoint,ClockWise> >
    : boost::mpl::int_<3>
{};

template <
    typename TPoint,
    bool ClockWise,
    std::size_t Index,
    std::size_t Dimension
>
struct indexed_access< Triangle<TPoint, ClockWise>, Index, Dimension>
{
    typedef typename geometry::coordinate_type<TPoint>::type coordinate_type;

    static inline
    coordinate_type get( Triangle<TPoint, ClockWise> const& triangle )
    {
        return geometry::get<Dimension>( triangle[Index] );
    }

    static inline
    void set( Triangle<TPoint, ClockWise> & triangle, coordinate_type const& value )
    {
        geometry::set<Dimension>( triangle[Index], value );
    }
};

template <
    typename TPoint,
    bool ClockWise
>
struct point_type< Triangle<TPoint, ClockWise> >
{
    typedef TPoint type;
};

}}} // namespace boost::geometry::traits

#endif // FLOE_GEOMETRY_GEOMETRIES_TRIANGLE_HPP

