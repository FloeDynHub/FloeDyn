#ifndef FLOE_GEOMETRY_GEOMETRIES_MULTI_CIRCLE_HPP
#define FLOE_GEOMETRY_GEOMETRIES_MULTI_CIRCLE_HPP

#include <memory>
#include <vector>

#include <boost/concept/assert.hpp>

#include "floe/geometry/core/tags.hpp"
#include "floe/geometry/geometries/concepts/multi_circle_concept.hpp"

namespace floe { namespace geometry
{

/*! MultiCircle, a collection of Circle
 */
template <
    typename TCircle,
    template <typename, typename> class TContainer = std::vector,
    template <typename> class TAllocator = std::allocator
>
struct MultiCircle
    : public TContainer<TCircle, TAllocator<TCircle> >
{
    BOOST_CONCEPT_ASSERT( (concept::Circle<TCircle>) );
};

}} // namespace floe::geometry

namespace boost { namespace geometry { namespace traits
{

//! Tag
template <
    typename TCircle,
    template <typename, typename> class TContainer,
    template <typename> class TAllocator
>
struct tag< floe::geometry::MultiCircle<TCircle, TContainer, TAllocator> >
{
    typedef multi_circle_tag type;
};

}}} // namespace boost::geometry::traits

#endif // FLOE_GEOMETRY_GEOMETRIES_MULTI_CIRCLE_HPP

