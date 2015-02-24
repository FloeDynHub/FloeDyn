/*!
 * \file integrate.hpp
 * \author Roland Denis
 */


#ifndef FLOE_INTEGRATION_INTEGRATE_HPP
#define FLOE_INTEGRATION_INTEGRATE_HPP

#include <utility>

#include "floe/integration/transform.hpp"

namespace floe { namespace integration
{

/*! Integration of a quadrature for an arbitrary shape
 *
 * \see transform
 */
template <
    typename TFunction,
    typename TGeometry,
    typename TStrategy
>
inline
typename floe::integration::detail::transform::return_type<TFunction,TGeometry>::type
    integrate( TFunction && function, TGeometry && geometry, TStrategy && strategy )
{
    return floe::integration::transform( std::forward<TFunction>(function), std::forward<TGeometry>(geometry), std::forward<TStrategy>(strategy) );
}

}} // namespace floe::integration

#endif // FLOE_INTEGRATION_INTEGRATE_HPP

