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
    using real_type = typename TPoint::value_type;

    Triangle( TPoint const& p1, TPoint const& p2, TPoint const& p3 )
        : std::template array<TPoint,3>{ { p1, p2, p3  }}
    {}


    public:
    inline real_type area()  { 
        TPoint& p1 {(*this)[0]}, p2 {(*this)[1]}, p3 {(*this)[2]};
        return abs((p3.x-p1.x)*(p2.y-p1.y)-(p2.x-p1.x)*(p3.y-p1.y)); 
    }

    inline bool is_fractured(TPoint crack_start,TPoint crack_stop)  { 
        TPoint& p1 {(*this)[0]}, p2 {(*this)[1]}, p3 {(*this)[2]};
        bool fract=false;
        real_type det= (p2.x-p1.x)*(crack_start.y-crack_stop.y)-(p2.y-p1.y)*(crack_start.x-crack_stop.x);
        real_type t1;
        if (abs(det)<0.01*((crack_start.y-crack_stop.y)*(crack_start.y-crack_stop.y)+(crack_start.x-crack_stop.x)*(crack_start.x-crack_stop.x))){
            t1=(crack_start.x-p1.x)*(crack_start.y-crack_stop.y)-(crack_start.x-crack_stop.x)*(crack_start.y-p1.y);
            if (t1/abs(det) >=0 || t1/abs(det)<=1){ fract=true;}
            t1=(p2.x-p1.x)*(crack_start.y-p1.y)-(crack_start.x-p1.x)*(p2.y-p1.y);
            if (t1/abs(det) >=0 || t1/abs(det)<=1){ fract=true;}
        }
        det= (p2.x-p3.x)*(crack_start.y-crack_stop.y)-(p2.y-p3.y)*(crack_start.x-crack_stop.x);
        if (abs(det)<0.01*((crack_start.y-crack_stop.y)*(crack_start.y-crack_stop.y)+(crack_start.x-crack_stop.x)*(crack_start.x-crack_stop.x))){
            t1=(crack_start.x-p3.x)*(crack_start.y-crack_stop.y)-(crack_start.x-crack_stop.x)*(crack_start.y-p3.y);
            if (t1/abs(det) >=0 || t1/abs(det)<=1){ fract=true;}
            t1=(p2.x-p3.x)*(crack_start.y-p3.y)-(crack_start.x-p3.x)*(p2.y-p3.y);
            if (t1/abs(det) >=0 || t1/abs(det)<=1){ fract=true;}
        }
        return fract; 
    }


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

