#ifndef FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_MULTI_SIMPLE_STATIC_POLYGON_CONCEPT_HPP
#define FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_MULTI_SIMPLE_STATIC_POLYGON_CONCEPT_HPP

#include <boost/concept_check.hpp>
#include <boost/range/concepts.hpp>
#include <boost/range/metafunctions.hpp>

#include "floe/geometry/geometries/concepts/simple_static_polygon_concept.hpp"

namespace boost { namespace geometry { namespace concepts
{

/*! multi simple static polygon concept
 *
 * See MultiPolygon.
 */
template <typename Geometry>
struct MultiSimpleStaticPolygon
{
    typedef typename boost::range_value<Geometry>::type polygon_type;

    BOOST_CONCEPT_ASSERT( (concepts::SimpleStaticPolygon<polygon_type>) );
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

    BOOST_CONCEPT_ASSERT( (concepts::ConstSimpleStaticPolygon<polygon_type>) );
    BOOST_CONCEPT_ASSERT( (boost::RandomAccessRangeConcept<Geometry>) );

public:

    BOOST_CONCEPT_USAGE(ConstMultiSimpleStaticPolygon)
    {
    }
};

}}} // namespace boost::geometry::concepts

namespace floe { namespace geometry { namespace concepts
{
    using boost::geometry::concepts::MultiSimpleStaticPolygon;
    using boost::geometry::concepts::ConstMultiSimpleStaticPolygon;

}}} // namespace floe::geometry::concepts

#endif // FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_MULTI_SIMPLE_STATIC_POLYGON_CONCEPT_HPP

