#ifndef FLOE_GEOMETRY_GEOMETRIES_POINT_HPP
#define FLOE_GEOMETRY_GEOMETRIES_POINT_HPP

#include <boost/mpl/int.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/coordinate_type.hpp>
#include <boost/geometry/core/coordinate_system.hpp>
#include <boost/geometry/core/coordinate_dimension.hpp>

namespace floe { namespace geometry {

/*! 2D Point type
 *
 * We could use point from boost::geometry. This is an exemple of what must be done to add our own type.
 * The coordinates are accessible using .x, .y, or with the [] operator, or using the values array.
 *
 * \tparam TCoordinate          Type of the coordinate
 * \tparam TCoordinateSystem    System of coordinate
 */
template <
    typename TCoordinate,
    typename TCoordinateSystem = boost::geometry::cs::cartesian
>
struct Point
{
    Point() {};
    
    Point( TCoordinate x, TCoordinate y ) : x{x}, y{y} {}   
   
    inline
    TCoordinate & operator[] (std::size_t i) {
        return values[i];
    }

    constexpr
    const TCoordinate & operator[] (std::size_t i) const {
        return values[i];
    }

    //! Access coordinate by name or index
    union {
        struct { TCoordinate x, y; };
        TCoordinate values[2];
    };

};

}} // namespace floe::geometry

// Type traits

// For easy use, locally
namespace {
    template <
        typename TCoordinate,
        typename TCoordinateSystem
    >
    using FloePoint = floe::geometry::Point<TCoordinate, TCoordinateSystem>;
}

// STL accessor
namespace std {

    template <
        size_t I,
        typename TCoordinate,
        typename TCoordinateSystem
    >
    constexpr
    TCoordinate & get( FloePoint<TCoordinate,TCoordinateSystem> & pt ) {
        return pt[I];
    }

    template <
        size_t I,
        typename TCoordinate,
        typename TCoordinateSystem
    >
    constexpr
    TCoordinate && get( FloePoint<TCoordinate,TCoordinateSystem> && pt ) {
        return pt[I];
    }

    template < 
        size_t I,
        typename TCoordinate,
        typename TCoordinateSystem
    >
    constexpr
    const TCoordinate & get( const FloePoint<TCoordinate,TCoordinateSystem> & pt ) {
        return pt[I];
    }

} // namespace std

// Boost geometry traits
namespace boost { namespace geometry { namespace traits
{


//! Tag
template <
    typename TCoordinate,
    typename TCoordinateSystem
>
struct tag< FloePoint< TCoordinate, TCoordinateSystem > >
{
    typedef point_tag type;
};

//! Coordinate Type
template <
    typename TCoordinate,
    typename TCoordinateSystem
>
struct coordinate_type< FloePoint< TCoordinate, TCoordinateSystem > >
{
    typedef TCoordinate type;
};

//! Coordinate system
template <
    typename TCoordinate,
    typename TCoordinateSystem
>
struct coordinate_system< FloePoint< TCoordinate, TCoordinateSystem > >
{
    typedef TCoordinateSystem type;
};


//! Dimension
template <
    typename TCoordinate,
    typename TCoordinateSystem
>
struct dimension< FloePoint< TCoordinate, TCoordinateSystem > >
    : boost::mpl::int_<2>
{};

//! Accessor
template <
    typename TCoordinate,
    typename TCoordinateSystem,
    std::size_t Dimension
>
struct access< FloePoint< TCoordinate, TCoordinateSystem >, Dimension >
{
    static inline TCoordinate get( const FloePoint<TCoordinate, TCoordinateSystem> & pt)
    {
        return pt[Dimension];
    }
        
    static inline void set(
        FloePoint<TCoordinate, TCoordinateSystem> & pt,
        TCoordinate const & value)
    {
        pt[Dimension] = value;
    }
};



}}} // namespace boost::geometry::traits


#endif // FLOE_GEOMETRY_GEOMETRIES_POINT_HPP
