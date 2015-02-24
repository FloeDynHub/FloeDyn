#ifndef FLOE_GEOMETRY_GEOMETRIES_POLYGON_HPP
#define FLOE_GEOMETRY_GEOMETRIES_POLYGON_HPP

#include <array>
#include "floe/geometry/geometries/ring.hpp"

namespace floe { namespace geometry {

/*! Polygon type
 *
 * /tparam TPoint       Point type
 * /tparam ClockWise    true for clockwise oriented polygon
 * /tparam Closed       true for closed polygon
 * /tparam TRing        Ring type
 */
template <
    typename TPoint,
    bool ClockWise  = true,
    bool Closed     = true,
    template< typename, bool, bool > class TRing  = Ring
>
class Polygon
{

public:

    // Member types
    using point_type = TPoint;
    using ring_type  = TRing<TPoint, ClockWise, Closed>;
    using inner_container_type = std::array<ring_type,0>; //! i.e. there is no inner ring

    //! Get the boundary
    inline       ring_type & boundary ()       { return m_boundary; }
    inline const ring_type & boundary () const { return m_boundary; }

private:
    ring_type m_boundary; //! Boundary of the polygon
};

}} // namespace floe::geometry

// Type traits

// For easy use, locally
namespace {
    template <
        typename TPoint,
        bool ClockWise,
        bool Closed,
        template < typename, bool, bool > class TRing
    >
    using FloePolygon = floe::geometry::Polygon<TPoint,ClockWise,Closed,TRing>;
} // local namespace

// Boost geometry traits
namespace boost { namespace geometry { namespace traits
{

//! Tag
template <
    typename TPoint,
    bool ClockWise,
    bool Closed,
    template < typename, bool, bool > class TRing
>
struct tag< FloePolygon<TPoint,ClockWise,Closed,TRing> >
{
    typedef polygon_tag type;
};


//! Ring const type
template <
    typename TPoint,
    bool ClockWise,
    bool Closed,
    template < typename, bool, bool > class TRing
>
struct ring_const_type< FloePolygon<TPoint,ClockWise,Closed,TRing> >
{
    typedef 
        typename FloePolygon<TPoint,ClockWise,Closed,TRing>::ring_type const&
        type;
};

//! Ring mutable type
template <
    typename TPoint,
    bool ClockWise,
    bool Closed,
    template < typename, bool, bool > class TRing
>
struct ring_mutable_type< FloePolygon<TPoint,ClockWise,Closed,TRing> >
{
    typedef 
        typename FloePolygon<TPoint,ClockWise,Closed,TRing>::ring_type &
        type;
};

//! Interior const type
template <
    typename TPoint,
    bool ClockWise,
    bool Closed,
    template < typename, bool, bool > class TRing
>
struct interior_const_type< FloePolygon<TPoint,ClockWise,Closed,TRing> >
{
    typedef 
        typename FloePolygon<TPoint,ClockWise,Closed,TRing>::inner_container_type
        type;
};

//! Interior mutable type
template <
    typename TPoint,
    bool ClockWise,
    bool Closed,
    template < typename, bool, bool > class TRing
>
struct interior_mutable_type< FloePolygon<TPoint,ClockWise,Closed,TRing> >
{
    typedef 
        typename FloePolygon<TPoint,ClockWise,Closed,TRing>::inner_container_type
        type;
};

//! Exterior ring
template <
    typename TPoint,
    bool ClockWise,
    bool Closed,
    template < typename, bool, bool > class TRing
>
struct exterior_ring< FloePolygon<TPoint,ClockWise,Closed,TRing> >
{
    typedef
        FloePolygon<TPoint,ClockWise,Closed,TRing>
        polygon_type;

    static inline
    typename polygon_type::ring_type& get( polygon_type& p ) {
        return p.boundary();
    }

    static inline
    typename polygon_type::ring_type const& get( polygon_type const& p ) {
        return p.boundary();
    }
};


//! Interior rings (none actually)
template <
    typename TPoint,
    bool ClockWise,
    bool Closed,
    template < typename, bool, bool > class TRing
>
struct interior_rings< FloePolygon<TPoint,ClockWise,Closed,TRing> >
{
    typedef
        FloePolygon<TPoint,ClockWise,Closed,TRing>
        polygon_type;

    static inline
    typename polygon_type::inner_container_type get( polygon_type& /* p */ ) {
        return {};
    }

    static inline
    typename polygon_type::inner_container_type get( polygon_type const& /* p */ ) {
        return {};
    }
};


}}} // namespace boost::geometry::traits



#endif // FLOE_GEOMETRY_GEOMETRIES_POLYGON_HPP
