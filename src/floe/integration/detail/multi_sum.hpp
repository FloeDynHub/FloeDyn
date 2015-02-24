/*!
 * \file multi_sum.hpp
 * \author Roland Denis
 *
 * Adapted from Boost.Geometry (aka GGL, Generic Geometry Library)
 */

#ifndef FLOE_INTEGRATION_DETAIL_MULTI_SUM_HPP
#define FLOE_INTEGRATION_DETAIL_MULTI_SUM_HPP

#include <boost/range.hpp>

namespace floe { namespace integration { namespace detail
{

/*! Sum of a strategy on a multi-geometry
 *
 * \see boost/geometry/multi/algorithms/detail/multi_sum.hpp
 */
struct multi_sum
{
    /*! Application
     *
     * \tparam TReturn      Return type
     * \tparam TPolicy      Type of the policy to be used on each geometry
     * \tparam TFunction    Type of the function to be integrate
     * \tparam TMultiGeometry   Type of the multi-geometry
     * \tparam TStrategy    Type of the strategy used for the integration
     */
    template <
        typename TReturn, 
        typename TPolicy, 
        typename TFunction,
        typename TMultiGeometry, 
        typename TStrategy
    >
    static inline TReturn apply(TFunction const& function, TMultiGeometry const& geometry, TStrategy const& strategy)
    {
        TReturn sum = TReturn();
        for (typename boost::range_iterator
                <
                    TMultiGeometry const
                >::type it = boost::begin(geometry);
            it != boost::end(geometry);
            ++it)
        {
            sum += TPolicy::apply(function, *it, strategy);
        }
        return sum;
    }
};

}}} // namespace boost::geometry::detail

#endif // FLOE_INTEGRATION_DETAIL_MULTI_SUM_HPP

