/*!
 * \file floe/collision/matlab/periodic_detector.hpp
 * \brief Collision detector (matlab version) and associated functions, taking space topology (periodicity) into account
 * \see MatlabDetector for more explanations.
 * \author Quentin Jouet
 */

#ifndef FLOE_COLLISION_MATLAB_PERIODIC_DETECTOR_DEF_HPP
#define FLOE_COLLISION_MATLAB_PERIODIC_DETECTOR_DEF_HPP

#include "floe/collision/matlab/periodic_detector.h"
#include "floe/collision/matlab/detector.hpp"


namespace floe { namespace collision { namespace matlab
{


//! Update detector
template <
    typename TFloeGroup,
    typename TSpaceTopology,
    typename TGhostFloe,
    typename TContact
>
void
PeriodicMatlabDetector<TFloeGroup, TSpaceTopology, TGhostFloe, TContact>::detect()
{
    // Number of floes
    const std::size_t N = base_class::get_nb_floes();
    const std::size_t Ng = base_class::m_prox_data.nb_ghosts();

    // Resize matrix
    base_class::m_prox_data.resize(N, N + Ng);

    // Detection
    base_class::detect_step1();
}


template <
    typename TFloeGroup,
    typename TSpaceTopology,
    typename TGhostFloe,
    typename TContact
>
typename PeriodicMatlabDetector<TFloeGroup, TSpaceTopology, TGhostFloe, TContact>::contact_type
PeriodicMatlabDetector<TFloeGroup, TSpaceTopology, TGhostFloe, TContact>::create_contact(
    std::size_t n1, std::size_t n2, point_type point1, point_type point2) const
{
    const std::size_t N { base_class::get_nb_floes() };
    if (n2 >= N){
        return { &base_class::get_floe(n1), &base_class::m_prox_data.get_ghost_floe(n2 - N), point1, point2 };
    } else if (n1 >= N) {
        return {&base_class::m_prox_data.get_ghost_floe(n1 - N), &base_class::get_floe(n2), point1, point2 };
    } else {
        return { &base_class::get_floe(n1), &base_class::get_floe(n2), point1, point2 };
    }
}


}}} // namespace floe::collision::matlab


#endif // FLOE_COLLISION_MATLAB_PERIODIC_DETECTOR_DEF_HPP