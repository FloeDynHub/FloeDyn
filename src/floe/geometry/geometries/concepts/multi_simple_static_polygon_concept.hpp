#ifndef FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_MULTI_SIMPLE_STATIC_POLYGON_CONCEPT_HPP
#define FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_MULTI_SIMPLE_STATIC_POLYGON_CONCEPT_HPP

#include <boost/concept_check.hpp>
#include <boost/range/concepts.hpp>
#include <boost/range/metafunctions.hpp>

#include "floe/geometry/geometries/concepts/simple_static_polygon_concept.hpp"

namespace boost { namespace geometry { namespace concept
{

/*! multi simple static polygon concept
 *
 * See MultiPolygon.
 */
template <typename Geometry>
struct MultiSimpleStaticPolygon
{
    typedef typename boost::range_value<Geometry>::type polygon_type;

    BOOST_CONCEPT_ASSERT( (concept::SimpleStaticPolygon<polygon_type>) );
    BOOST_CONCEPT_ASSERT( (boost::RandomAccessRangeConcept<Geometry>) );

public:

    BOOST_CONCEPT_USAGE(MultiSimpleStaticPolygon)
    {
    }
};

/*! multi simple static polygon concept (const version)
 *
 * See ConstMultiPolygon.
 */
template <typename Geometry>
struct ConstMultiSimpleStaticPolygon
{
    typedef typename boost::range_value<Geometry>::type polygon_type;

    BOOST_CONCEPT_ASSERT( (concept::ConstSimpleStaticPolygon<polygon_type>) );
    BOOST_CONCEPT_ASSERT( (boost::RandomAccessRangeConcept<Geometry>) );

public:

    BOOST_CONCEPT_USAGE(ConstMultiSimpleStaticPolygon)
    {
    }
};

}}} // namespace boost::geometry::concept

namespace floe { namespace geometry { namespace concept
{
    using boost::geometry::concept::MultiSimpleStaticPolygon;
    using boost::geometry::concept::ConstMultiSimpleStaticPolygon;

}}} // namespace floe::geometry::concept

#endif // FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_MULTI_SIMPLE_STATIC_POLYGON_CONCEPT_HPP

