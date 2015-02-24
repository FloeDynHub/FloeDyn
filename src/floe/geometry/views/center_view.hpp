#ifndef FLOE_GEOMETRY_VIEWS_CENTER_VIEW_HPP
#define FLOE_GEOMETRY_VIEWS_CENTER_VIEW_HPP

//! Copy of boost/geometry/extensions/nsphere/views/center_view.hpp

#include <boost/geometry/core/point_type.hpp>

namespace boost { namespace geometry
{

template < typename TCircle >
struct center_view
{
    typedef typename geometry::point_type<TCircle>::type      point_type;
    typedef typename geometry::coordinate_type<TCircle>::type coordinate_type;

    explicit center_view( TCircle & circle )
        : m_circle(circle)
    {}

    template <std::size_t I> coordinate_type get() const { return geometry::get<I>(m_circle); }
    template <std::size_t I> void set( coordinate_type const& value ) { geometry::set<I>(m_circle, value); }
    
private:
    TCircle & m_circle;
};

// Type traits
namespace traits
{

template < typename TCircle >
struct tag< center_view<TCircle> >
{
    typedef point_tag type;
};

template < typename TCircle >
struct coordinate_type< center_view<TCircle> >
{
    typedef 
        typename geometry::coordinate_type< typename geometry::point_type<TCircle>::type >::type type;
};

template < typename TCircle >
struct coordinate_system< center_view<TCircle> >
{
    typedef
        typename geometry::coordinate_system< typename geometry::point_type<TCircle>::type >::type type;
};

template < typename TCircle >
struct dimension< center_view<TCircle> >
    : geometry::dimension< typename geometry::point_type<TCircle>::type >
{};

template <
    typename    TCircle,
    std::size_t Dimension
>
struct access< center_view<TCircle>, Dimension >
{
    typedef
        typename geometry::coordinate_type< typename geometry::point_type<TCircle>::type >::type
        coordinate_type;

    static inline
    coordinate_type get( center_view<TCircle> const& center )
    {
        return center.template get<Dimension>();
    }

    static inline
    void set( center_view<TCircle> & center, coordinate_type const& value )
    {
        center.template set<Dimension>(value);
    }
};


} // namespace traits


}} // namespace boost::geometry

namespace floe { namespace geometry
{

using boost::geometry::center_view;

}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_VIEWS_CENTER_VIEW_HPP
