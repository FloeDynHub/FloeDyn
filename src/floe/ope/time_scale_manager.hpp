/*!
 * \file ope/time_scale_manager.hpp
 * \brief Time Scale Manager
 * \author Quentin Jouet
 */

#ifndef OPE_TIME_SCALE_MANAGER_HPP
#define OPE_TIME_SCALE_MANAGER_HPP

#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/multi_point.hpp"
#include "floe/geometry/frame/frame_transformers.hpp"

#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace floe { namespace ope
{

/*! TimeScaleManager
 *
 * Operator for time step determination
 * i.e. finding maximal time step that ensures no floes interpenetration
 *
 */

template<typename TDetector>
class TimeScaleManager
{

public:
    using value_type = typename TDetector::value_type;
    using point_type = typename TDetector::point_type;
    using floe_type = typename TDetector::floe_type;
    using optim_type = typename TDetector::optim_type;
    using frame_type = typename floe_type::frame_type;
    using multi_point_type = floe::geometry::MultiPoint<point_type>;
    using floe_interface_type = typename TDetector::floe_interface_type;
    using optim_interface_type = typename TDetector::optim_interface_type;

    /*!
     * Returns time step taking all floes into account, base on detector informations
     */
    template <typename TDomain>
    value_type delta_t_secu(TDomain* domain, TDetector* detector);

private:

    /*!
     * Returns maximal delta_t for 2 floes beeing close
     * Corresponds to gestion_temps() in Matlab code
     */
    value_type delta_t_secu(
        value_type dist_secu,
        value_type dist_opt,
        const floe_type& floe1,
        const floe_interface_type& floe2,
        const optim_type& optim1,
        const optim_interface_type& optim2,
        short I,
        value_type dt_default
    );

    /*!
     * Returns maximal delta_t for 2 floes not beeing close
     * Corresponds to gestion_temps_fast() in Matlab code
     */
    value_type delta_t_secu_fast(
        value_type dist_secu,
        const floe_type& floe1,
        const floe_interface_type& floe2,
        const optim_type& optim1,
        const optim_interface_type& optim2
    );

};

template <typename TDetector>
template <typename TDomain>
typename TimeScaleManager<TDetector>::value_type
TimeScaleManager<TDetector>::delta_t_secu(TDomain* domain, TDetector* detector)
{
    value_type dt_default = domain->default_time_step();
    
    if (detector->interpenetration())
    {
        domain->set_time_step(domain->time_step() / 2);
        return domain->time_step();
    }

    detector->clean_dist_opt();
    auto const& dist_secu = detector->m_dist_secu;
    auto& dist_opt = detector->m_dist_opt;
    auto& indics = detector->m_indic;
    auto& floes = detector->m_floes;
    auto& optims = detector->m_optims;
    value_type global_min_dt = dt_default;
    
    #pragma omp parallel for
    for (std::size_t i = 0; i < dist_secu.size1(); ++i)
    {
        for ( std::size_t j = i+ 1; j != dist_secu.size2(); ++j )
        {
            value_type delta_t;
            if (indics(i,j) == 0)
            {   
                delta_t = delta_t_secu_fast(dist_secu(i,j),
                    *floes[i], detector->get_floe(j), *optims[i], detector->get_optim(j));
            }
            else
            {
                delta_t = delta_t_secu(dist_secu(i,j), dist_opt(i,j),
                         *floes[i], detector->get_floe(j), *optims[i], detector->get_optim(j), indics(i,j), dt_default);  
            }

            global_min_dt = std::min(
                global_min_dt,
                delta_t
            );
        }
    }


    domain->set_time_step(global_min_dt);
    return global_min_dt;
}

template <typename TDetector>
typename TimeScaleManager<TDetector>::value_type
TimeScaleManager<TDetector>::delta_t_secu(
    value_type dist_secu,
    value_type dist_opt,
    const floe_type& floe1,
    const floe_interface_type& floe2,
    const optim_type& optim1,
    const optim_interface_type& optim2,
    short I,
    value_type dt_default
){
    const point_type& Vg1 = floe1.state().speed, Vg2 = floe2.state().speed;
    const value_type& Vt1 = floe1.state().rot, Vt2 = floe2.state().rot;
    const point_type& C1 = optim1.global_disk().center, C2 = optim2.global_disk().center;
    const value_type& R1 = optim1.global_disk().radius, R2 = optim2.global_disk().radius;
    const point_type& G1 = floe1.state().pos, G2 = floe2.state().pos;
    const value_type& tau1 = optim1.tau(), tau2 = optim2.tau();
    const value_type& dc1 = optim1.cdist(), dc2 = optim2.cdist();
    value_type lambda = std::min(dc1, dc2) / 20;

    value_type d = std::max(dist_opt, dist_secu);
    lambda = std::min(lambda, d / 20);

    // Calcul du deplacement d un point par rapport aux reperes en t+dt.
    // repere a l'instant t+dt_default :
    frame_type mark1{G1 + dt_default * (Vg1 - Vg2), dt_default * Vt1};
    frame_type mark2{G2 + dt_default * (Vg2 - Vg1), dt_default * Vt2};
    // fin calcul repere

    // calcul ceinture de points: (pour la rotation)
    auto D1 = (distance(C1, G1) + R1 - tau1);
    auto D2 = (distance(C2,G2) + R2 - tau2);
    multi_point_type Belt_P1, Belt_P2;
    for (int i=0; i<50; ++i)
    {
        auto angle = 2 * i *  M_PI / 50;
        Belt_P1.push_back(point_type{D1 * cos(angle), D1 * sin(angle)});
        Belt_P2.push_back(point_type{D2 * cos(angle), D2 * sin(angle)});
    }

    multi_point_type Belt_P1_be, Belt_P2_be, Belt_P1_af, Belt_P2_af;

    using namespace geometry::frame; // import transformer, itransformer
    /*
    // version matlab raccourcie
    geometry::transform( Belt_P1, Belt_P1_af, transformer( frame_type{C1 - G1, 0} ));
    geometry::transform( Belt_P2, Belt_P2_af, transformer( frame_type{C2 - G2, 0} ));
    geometry::transform( Belt_P1_af, Belt_P1_af, transformer( mark1, mark2 ));
    geometry::transform( Belt_P2_af, Belt_P2_af, transformer( mark2, mark1 ));

    value_type dist1 = 0, dist2 = 0;
    for (std::size_t i = 0; i != Belt_P1.size(); ++i)
        dist1 = std::max(dist1, distance(Belt_P1[i] + C1, Belt_P1_af[i] + G2));
    for (std::size_t i = 0; i != Belt_P2.size(); ++i)
        dist2 = std::max(dist2, distance(Belt_P2[i] + C2, Belt_P2_af[i] + G1));
    // version raccourcie
    */

    // Version 2
    geometry::transform( Belt_P1, Belt_P1_be, transformer( frame_type{G1, 0} ));
    geometry::transform( Belt_P2, Belt_P2_be, transformer( frame_type{G2, 0} ));
    geometry::transform( Belt_P1, Belt_P1_af, transformer( mark1 ));
    geometry::transform( Belt_P2, Belt_P2_af, transformer( mark2 ));
    
    value_type dist1 = 0, dist2 = 0;
    for (std::size_t i = 0; i != Belt_P1.size(); ++i)
        dist1 = std::max(dist1, distance(Belt_P1_be[i], Belt_P1_af[i]));
    for (std::size_t i = 0; i != Belt_P2.size(); ++i)
        dist2 = std::max(dist2, distance(Belt_P2_be[i], Belt_P2_af[i]));
    // END Version 2
    

    // %%%%%%%%%% Cas obj1 %%%%%%%%%%
    value_type calc;
    if (dist1 < 1e-12)
        calc = dt_default;
    else
        calc = ((d - lambda) / 2) * (dt_default / dist1);
    value_type dt12 = std::min(dt_default, calc);
    // %%%%%%%%%% Fin %%%%%%%%%%
    
    // %%%%%%%%%% Cas obj2 %%%%%%%%%%
    if (dist2 < 1e-12)
        calc = dt_default;
    else
        calc = ((d - lambda) / 2) * (dt_default / dist2);
    value_type dt21 = std::min(dt_default, calc);
    // %%%%%%%%%% Fin %%%%%%%%%%

    assert( std::min(dt12, dt21) > 0 );

    return std::min(dt12, dt21);
}

template <typename TDetector>
typename TimeScaleManager<TDetector>::value_type
TimeScaleManager<TDetector>::delta_t_secu_fast(
        value_type dist_secu,
        const floe_type& floe1,
        const floe_interface_type& floe2,
        const optim_type& optim1,
        const optim_interface_type& optim2
){

    point_type Vg1 = floe1.state().speed, Vg2 = floe2.state().speed;
    point_type C1 = optim1.global_disk().center, C2 = optim2.global_disk().center;
    value_type dc1 = optim1.cdist(), dc2 = optim2.cdist();

    // Calcul de la marge lambda
    value_type lambda = std::min(dc1, dc2) / 20;
    
    // Axe reliant chaque couple de floe
    point_type Axe = (C2 - C1) / distance(C1, C2);

    // Vitesse relative projetÃ©e sur l'axe
    value_type VRel = geometry::dot_product(Vg2 - Vg1, Axe);

    // Calcul du dt max
    value_type delta_t;
    
    if (VRel < 0)
    {
        // Collision possible
        delta_t = - ( dist_secu - lambda ) / VRel;
    } else
    {
        // Collision impossible
        delta_t = std::numeric_limits<value_type>::max();
    }

    return delta_t;

}


}} // namespace floe::ope


#endif // OPE_TIME_SCALE_MANAGER_HPP