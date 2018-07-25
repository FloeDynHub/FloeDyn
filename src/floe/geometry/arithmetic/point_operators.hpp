#ifndef FLOE_GEOMETRY_ARITHMETIC_POINT_OPERATORS
#define FLOE_GEOMETRY_ARITHMETIC_POINT_OPERATORS

#include <cmath>
#include <cstddef>     // std::size_t
#include <type_traits> // std::enable_if
#include <ostream>
#include "floe/geometry/geometries/point.hpp"

namespace floe { namespace geometry
{

// Oulà ...
#include "floe/arithmetic/container_operators.hpp"

//namespace floe { namespace arithmetic { namespace container_operators
//{
namespace container_operators
{

namespace fg = ::floe::geometry;

//! Activate operators for point type
template <
    typename TCoordinate,
    typename TCoordinateSystem
>
struct is_activated< fg::Point<TCoordinate, TCoordinateSystem> >
    : yes {};

// Traits

//! Traits for size
template <
    typename TCoordinate,
    typename TCoordinateSystem
>
struct size_impl< fg::Point<TCoordinate, TCoordinateSystem> >
{
    static constexpr size_t size()
    {
        return 2;
    }
};

//! Traits for constructor
template <
    typename TCoordinate,
    typename TCoordinateSystem
>
struct cstr_impl< fg::Point<TCoordinate, TCoordinateSystem> >
{
    template< typename V1, typename V2 >
    static constexpr
    fg::Point<TCoordinate, TCoordinateSystem>
        from_values( V1 x, V2 y )
    {
        return {x, y};
    }
};

//! Traits for value type
template <
    typename TCoordinate,
    typename TCoordinateSystem
>
struct valueType< fg::Point<TCoordinate, TCoordinateSystem> >
{
    typedef TCoordinate type;
};

//}}} // namespace floe::arithmetic::container_operators

}}} // namespace floe::geometry::container_operators


#endif // FLOE_GEOMETRY_ARITHMETIC_POINT_OPERATORS
