#ifndef FLOE_GEOMETRY_GEOMETRIES_CONCEPT_CHECK_HPP
#define FLOE_GEOMETRY_GEOMETRIES_CONCEPT_CHECK_HPP

#include <boost/geometry/geometries/concepts/check.hpp>

#include "floe/geometry/core/tag.hpp"
#include "floe/geometry/core/tags.hpp"

#include "floe/geometry/geometries/concepts/circle_concept.hpp"
#include "floe/geometry/geometries/concepts/static_ring_concept.hpp"
#include "floe/geometry/geometries/concepts/static_polygon_concept.hpp"
#include "floe/geometry/geometries/concepts/simple_static_polygon_concept.hpp"
#include "floe/geometry/geometries/concepts/triangle_concept.hpp"
#include "floe/geometry/geometries/concepts/mesh_concept.hpp"
#include "floe/geometry/geometries/concepts/multi_simple_static_polygon_concept.hpp"
#include "floe/geometry/geometries/concepts/multi_point_concept.hpp"
#include "floe/geometry/geometries/concepts/multi_circle_concept.hpp"

namespace boost { namespace geometry {

namespace dispatch {


/*
template <
    typename Geometry,
    bool IsConst
>
struct check<Geometry, circle_tag, IsConst>
{};
*/

template <typename Geometry>
struct check< Geometry, circle_tag, true >
    : detail::concept_check::check< concept::ConstCircle< Geometry > >
{};

template < typename Geometry >
struct check< Geometry, circle_tag, false >
    : detail::concept_check::check< concept::Circle< Geometry > >
{};

template < typename Geometry >
struct check< Geometry, multi_circle_tag, true >
    : detail::concept_check::check< concept::ConstMultiCircle<Geometry> >
{};

template < typename Geometry >
struct check< Geometry, multi_circle_tag, false >
    : detail::concept_check::check< concept::MultiCircle<Geometry> >
{};

template < typename Geometry >
struct check< Geometry, static_polygon_tag, true >
    : detail::concept_check::check< concept::ConstStaticPolygon< Geometry > >
{};

template < typename Geometry >
struct check< Geometry, static_polygon_tag, false >
    : detail::concept_check::check< concept::StaticPolygon<Geometry> >
{};

template < typename Geometry >
struct check< Geometry, static_ring_tag, true >
    : detail::concept_check::check< concept::ConstStaticRing<Geometry> >
{};

template < typename Geometry >
struct check< Geometry, static_ring_tag, false >
    : detail::concept_check::check< concept::StaticRing<Geometry> >
{};


template < typename Geometry >
struct check< Geometry, triangle_tag, true >
    : detail::concept_check::check< concept::ConstTriangle<Geometry> >
{};

template < typename Geometry >
struct check< Geometry, triangle_tag, false >
    : detail::concept_check::check< concept::Triangle<Geometry> >
{};


template < typename Geometry >
struct check< Geometry, simple_static_polygon_tag, true >
    : detail::concept_check::check< concept::ConstSimpleStaticPolygon<Geometry> >
{};

template < typename Geometry >
struct check< Geometry, simple_static_polygon_tag, false >
    : detail::concept_check::check< concept::SimpleStaticPolygon<Geometry> >
{};


template < typename Geometry >
struct check< Geometry, multi_simple_static_polygon_tag, true >
    : detail::concept_check::check< concept::ConstMultiSimpleStaticPolygon<Geometry> >
{};

template < typename Geometry >
struct check< Geometry, multi_simple_static_polygon_tag, false >
    : detail::concept_check::check< concept::MultiSimpleStaticPolygon<Geometry> >
{};


template < typename Geometry >
struct check< Geometry, mesh_tag, true >
    : detail::concept_check::check< concept::ConstMesh<Geometry> >
{};

template < typename Geometry >
struct check< Geometry, mesh_tag, false >
    : detail::concept_check::check< concept::Mesh<Geometry> >
{};

} // namespace dispatch

}} // namespace boost::geometry

namespace floe { namespace geometry { namespace concept
{
    using boost::geometry::concept::check;
}}} // namespace floe::geometry::concept

#endif // FLOE_GEOMETRY_GEOMETRIES_CONCEPT_CHECK_HPP
