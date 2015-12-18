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
    typename TFloe,
    typename TSpaceTopology,
    typename TGhostFloe,
    typename TContact
>
void
PeriodicMatlabDetector<TFloe, TSpaceTopology, TGhostFloe, TContact>::detect()
{
    // Number of floes
    const std::size_t N = base_class::m_floes.size();
    const std::size_t Ng = m_ghost_floes.size();

    // Resize matrix
    base_class::m_indic.resize(N, N + Ng);
    base_class::m_dist_opt = ublas::scalar_matrix<value_type>(N, N + Ng, 0);
    base_class::m_dist_secu = ublas::scalar_matrix<value_type>(N, N + Ng, 0);

    // Detection
    base_class::detect_step1();
}


template <
    typename TFloe,
    typename TSpaceTopology,
    typename TGhostFloe,
    typename TContact
>
typename PeriodicMatlabDetector<TFloe, TSpaceTopology, TGhostFloe, TContact>::floe_interface_type const&
PeriodicMatlabDetector<TFloe, TSpaceTopology, TGhostFloe, TContact>::get_floe(std::size_t n) const
{
    const std::size_t N { base_class::m_floes.size() };
    if (n < N)
        return *(base_class::m_floes[n]);
    else
        return m_ghost_floes[n - N];
}


template <
    typename TFloe,
    typename TSpaceTopology,
    typename TGhostFloe,
    typename TContact
>
typename PeriodicMatlabDetector<TFloe, TSpaceTopology, TGhostFloe, TContact>::optim_interface_type const&
PeriodicMatlabDetector<TFloe, TSpaceTopology, TGhostFloe, TContact>::get_optim(std::size_t n) const
{
    const std::size_t N { base_class::m_floes.size() };
    if (n < N)
        return *(base_class::m_optims[n]);
    else
        return m_ghost_optims[n - N];
}


template <
    typename TFloe,
    typename TSpaceTopology,
    typename TGhostFloe,
    typename TContact
>
void
PeriodicMatlabDetector<TFloe, TSpaceTopology, TGhostFloe, TContact>::set_dist_opt(
    std::size_t n1, std::size_t n2, value_type val)
{
    const std::size_t N { base_class::m_floes.size() };
    if (n2 >= N){
        base_class::m_dist_opt(n1, n2) = val;
    } else if (n1 >= N) {
        base_class::m_dist_opt(n2, n1) = val;
    } else {
        base_class::m_dist_opt(n1, n2) = base_class::m_dist_opt(n2, n1) = val;
    }
}


template <
    typename TFloe,
    typename TSpaceTopology,
    typename TGhostFloe,
    typename TContact
>
typename PeriodicMatlabDetector<TFloe, TSpaceTopology, TGhostFloe, TContact>::contact_type
PeriodicMatlabDetector<TFloe, TSpaceTopology, TGhostFloe, TContact>::create_contact(
    std::size_t n1, std::size_t n2, point_type point1, point_type point2) const
{
    const std::size_t N { base_class::m_floes.size() };
    if (n2 >= N){
        return {base_class::m_floes[n1], &m_ghost_floes[n2 - N], point1, point2 };
    } else if (n1 >= N) {
        return {&m_ghost_floes[n1 - N], base_class::m_floes[n2], point1, point2 };
    } else {
        return { base_class::m_floes[n1], base_class::m_floes[n2], point1, point2 };
    }
}


}}} // namespace floe::collision::matlab


#endif // FLOE_COLLISION_MATLAB_PERIODIC_DETECTOR_DEF_HPP