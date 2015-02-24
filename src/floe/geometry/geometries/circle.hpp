#ifndef FLOE_GEOMETRY_GEOMETRIES_CIRCLE_HPP
#define FLOE_GEOMETRY_GEOMETRIES_CIRCLE_HPP

#include <boost/concept/assert.hpp>
#include <boost/geometry/geometries/concepts/point_concept.hpp>

namespace floe { namespace geometry {

/*! Circle type
 *
 * \tparam TPoint   Center point type
 */
template <
    typename TPoint,
    typename TRadius = typename boost::geometry::traits::coordinate_type<TPoint>::type
>
class Circle
{
    BOOST_CONCEPT_ASSERT( (boost::geometry::concept::Point<TPoint>) );

public:

    // Member types
    using point_type = TPoint;
    using coordinate_type = typename boost::geometry::traits::coordinate_type<TPoint>::type;
    using radius_type = TRadius;

    Circle() : center{0,0}, radius{0} {}
    explicit Circle(radius_type r) : center{{0,0}}, radius{r} {}
    Circle(point_type const& c, radius_type r) : center{c}, radius{r} {}

    point_type      center;
    radius_type     radius;

private:

}; // class Circle

}} // namespace floe::geometry


// Type traits

// For easy use, locally
namespace {
    template <
        typename TPoint,
        typename TRadius
    >
    using FloeCircle = floe::geometry::Circle<TPoint, TRadius>;
} // local namespace

// Boost geometry traits
namespace boost { namespace geometry { namespace traits {

//! Tag
template < typename TPoint, typename TRadius >
struct tag< FloeCircle<TPoint, TRadius> >
{
    typedef circle_tag type;
};

//! Coordinate Type
template < typename TPoint, typename TRadius >
struct coordinate_type< FloeCircle<TPoint, TRadius> >
{
    typedef typename FloeCircle<TPoint, TRadius>::coordinate_type type;
};

//! Coordinate System
template < typename TPoint, typename TRadius >
struct coordinate_system< FloeCircle<TPoint, TRadius> >
{
    typedef typename coordinate_system<TPoint>::type type;
};

//! Center Type
template < typename TPoint, typename TRadius >
struct point_type< FloeCircle<TPoint, TRadius> >
{
    typedef TPoint type;
};

//! Radius type
template < typename TPoint, typename TRadius >
struct radius_type< FloeCircle<TPoint, TRadius> >
{
    typedef TRadius type;
};

//! Center accessor
template <
    typename TPoint,
    typename TCoordinate,
    std::size_t Dimension
>
struct access< FloeCircle<TPoint, TCoordinate>, Dimension>
{
    static inline
    TCoordinate get( FloeCircle<TPoint, TCoordinate> const& circle )
    {
        return geometry::get<Dimension>(circle.center);
    }

    static inline
    void set( FloeCircle<TPoint, TCoordinate> & circle, TCoordinate const& value )
    {
        geometry::set<Dimension>(circle.center, value);
    }
};

//! Radius accessor
template <
    typename TPoint,
    typename TRadius
>
struct radius_access< FloeCircle<TPoint, TRadius>, TRadius>
{
    static inline
    TRadius get( FloeCircle<TPoint, TRadius> const& circle )
    {
        return circle.radius;
    }

    static inline
    void set( FloeCircle<TPoint, TRadius> & circle, TRadius const& value )
    {
        circle.radius = value;
    }
};


}}} // namespace boost::geometry::traits

#endif // FLOE_GEOMETRY_GEOMETRIES_CIRCLE_HPP
