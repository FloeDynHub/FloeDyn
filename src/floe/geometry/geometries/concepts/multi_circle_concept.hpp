#ifndef FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_MULTI_CIRCLE_HPP
#define FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_MULTI_CIRCLE_HPP

#include <boost/concept_check.hpp>
#include <boost/range/concepts.hpp>
#include <boost/range/metafunctions.hpp>

#include "floe/geometry/geometries/concepts/circle_concept.hpp"

namespace boost { namespace geometry { namespace concept
{

/*! multi circle concept
 *
 * See other multi concepts.
 */
template <typename Geometry>
struct MultiCircle
{
    typedef typename boost::range_value<Geometry>::type circle_type;

    BOOST_CONCEPT_ASSERT( (concept::Circle<circle_type>) );
    BOOST_CONCEPT_ASSERT( (boost::RandomAccessRangeConcept<Geometry>) );

public:
    
    BOOST_CONCEPT_USAGE(MultiCircle)
    {
    }
};


/*! multi circle concept (const version)
 *
 * See other const multi concepts.
 */
template <typename Geometry>
struct ConstMultiCircle
{
    typedef typename boost::range_value<Geometry>::type circle_type;

    BOOST_CONCEPT_ASSERT( (concept::ConstCircle<circle_type>) );
    BOOST_CONCEPT_ASSERT( (boost::RandomAccessRangeConcept<Geometry>) );

public:
    
    BOOST_CONCEPT_USAGE(ConstMultiCircle)
    {
    }
};

}}} // namespace boost::geometry::concepts

namespace floe { namespace geometry { namespace concept
{
    using boost::geometry::concepts::MultiCircle;
    using boost::geometry::concepts::ConstMultiCircle;
}}} // namespace floe::geometry::concept

#endif // FLOE_GEOMETRY_GEOMETRIES_CONCEPTS_MULTI_CIRCLE_HPP

