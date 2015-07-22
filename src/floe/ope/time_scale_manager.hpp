/*!
 * \file ope/time_scale_manager.hpp
 * \brief time_scale manager
 * \author Quentin Jouet
 */

#ifndef OPE_TIME_SCALE_MANAGER_HPP
#define OPE_TIME_SCALE_MANAGER_HPP

#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/multi_point.hpp"
#include "floe/geometry/frame/frame_transformers.hpp"

#include <math.h>

#include <iostream> // DEBUG

namespace floe { namespace ope
{

/*! TimeScaleManager
 *
 * Operator for time step determination
 *
 */


template <typename TDomain, typename TDetector>
class TimeScaleManager
{

public:

    // using floe_group_type = TFloeGroup;
    using value_type = typename TDetector::value_type;
    using point_type = typename TDetector::point_type;
    using floe_type = typename TDetector::floe_type;
    using optim_type = typename TDetector::optim_type;
    using frame_type = typename floe_type::frame_type;
    using multi_point_type = floe::geometry::MultiPoint<point_type>;
    using floe_interface_type = typename TDetector::floe_interface_type;
    using optim_interface_type = typename TDetector::optim_interface_type;

    // Default constructor
    TimeScaleManager() : m_domain{nullptr} {} //, m_detector{nullptr} {}

    // TimeScaleManager(TDomain& domain, TDetector& detector) : m_domain{&domain}, m_detector{&detector} {}

    value_type delta_t_secu(TDomain* domain, TDetector* detector);

private:

    TDomain* m_domain;
    // TDetector* m_detector;

    //! Corresponds to gestion_temps() in Matlab code
    value_type delta_t_secu(
        std::size_t I,
        value_type dist_secu,
        value_type dist_opt,
        const floe_type& floe1,
        const floe_interface_type& floe2,
        const optim_type& optim1,
        const optim_interface_type& optim2
        // + ghost floe cell translation
    );

    //! Corresponds to gestion_temps_fast() in Matlab code
    value_type delta_t_secu_fast(
        value_type dist_secu,
        const floe_type& floe1,
        const floe_interface_type& floe2,
        const optim_type& optim1,
        const optim_interface_type& optim2
        // + ghost floe cell translation
    );

};

template <typename TDomain, typename TDetector>
typename TimeScaleManager<TDomain, TDetector>::value_type
TimeScaleManager<TDomain, TDetector>::delta_t_secu(TDomain* domain, TDetector* m_detector)
{
    m_domain = domain;
    auto& dist_secu = m_detector->get_dist_secu();
    auto& dist_opt = m_detector->get_dist_opt();
    auto& indics = m_detector->get_indic();
    auto& floes = m_detector->get_floes();
    auto& optims = m_detector->get_optims();

    value_type global_min_dt = m_domain->default_time_step();
    // int nb{0}, nb_f{0};

    for (std::size_t i = 0; i!= dist_secu.size1(); ++i)
    {
        for ( std::size_t j = i+ 1; j != dist_secu.size2(); ++j )
        {
            if (indics(i,j) == 0)
            {
                global_min_dt = std::min(
                    global_min_dt,
                    delta_t_secu_fast(dist_secu(i,j),
                         *floes[i], m_detector->get_floe(j), *optims[i], m_detector->get_optim(j))
                ); //nb_f++;
            }
            else
            {
                global_min_dt = std::min(
                    global_min_dt,
                    delta_t_secu(indics(i,j), dist_secu(i,j), dist_opt(i,j),
                         *floes[i], m_detector->get_floe(j), *optims[i], m_detector->get_optim(j))
                ); //nb++;
            }
        }
    }

    // std::cout << "nb call dt_secu: " << nb << std::endl; // DEBUG

    // return std::max(global_min_dt, 1e-3);
    return global_min_dt;
}


template <typename TDomain, typename TDetector>
typename TimeScaleManager<TDomain, TDetector>::value_type
TimeScaleManager<TDomain, TDetector>::delta_t_secu(
    std::size_t I,
    value_type dist_secu,
    value_type dist_opt,
    const floe_type& floe1,
    const floe_interface_type& floe2,
    const optim_type& optim1,
    const optim_interface_type& optim2
){
    const value_type& dt_defaut = m_domain->default_time_step();
    const point_type& Vg1 = floe1.state().speed, Vg2 = floe2.state().speed;
    const value_type& Vt1 = floe1.state().rot, Vt2 = floe2.state().rot;
    const point_type& C1 = optim1.global_disk().center, C2 = optim2.global_disk().center;
    const value_type& R1 = optim1.global_disk().radius, R2 = optim2.global_disk().radius;
    const point_type& G1 = floe1.state().pos, G2 = floe2.state().pos;
    const value_type& tau1 = optim1.tau(), tau2 = optim2.tau();
    const value_type& dc1 = optim1.cdist(), dc2 = optim2.cdist();
    value_type lambda = std::min(dc1, dc2) / 20;

    value_type d;
    if (dist_opt != 0) { d = dist_opt; } else { d = dist_secu; }
    lambda = std::min(lambda, d / 20); // TODO éclaircir avec Mathias (lambda > d dans certains cas)

    // if I == 0  le pas de temps est calcule dans gestion_temps_fast.m
    // On est dans le cas I=1 :
    // Calcul du deplacement d un point par rapport aux reperes en t+dt.

    // repere a l'instant t+dt_defaut :
    frame_type mark1{G1 + dt_defaut * Vg1, dt_defaut * Vt1};
    frame_type mark2{G2 + dt_defaut * Vg2, dt_defaut * Vt2};
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
    /* Equivalent Matlab (trop lent)
    //passage dans repabs(t), puis rep1(t) et enfin repabs(t+dt):
    geometry::transform( Belt_P1, Belt_P1_be, transformer( frame_type{C1, 0} ));
    geometry::transform( Belt_P1_be, Belt_P1, itransformer( frame_type{G1, 0} ));
    geometry::transform( Belt_P2, Belt_P2_be, transformer( frame_type{C2, 0} ));
    geometry::transform( Belt_P2_be, Belt_P2, itransformer( frame_type{G2, 0} ));
    // passage dans rep2(t+dt) puis repabs(t):
    geometry::transform( Belt_P1, Belt_P1, transformer( mark1, mark2 ));
    geometry::transform( Belt_P2, Belt_P2, transformer( mark2, mark1 ));
    geometry::transform( Belt_P1, Belt_P1_af, transformer( frame_type{G2, 0} ));
    geometry::transform( Belt_P2, Belt_P2_af, transformer( frame_type{G1, 0} ));

    // calcul distance parcourue:
    value_type dist1 = 0, dist2 = 0;
    for (std::size_t i = 0; i != Belt_P1.size(); ++i)
        dist1 = std::max(dist1, distance(Belt_P1_be[i], Belt_P1_af[i]));
    for (std::size_t i = 0; i != Belt_P2.size(); ++i)
        dist2 = std::max(dist2, distance(Belt_P2_be[i], Belt_P2_af[i]));
    */
    /*
    // version raccourcie
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
    
    // version raccourcie + modifiée car transform(A, B, strategy) anormalement lent
    multi_point_type Belt_P1_copy{Belt_P1}, Belt_P2_copy{Belt_P2}; // waiting for transform optimization !
    geometry::transform( Belt_P1, Belt_P1, transformer( frame_type{C1 - G1, 0} ));
    geometry::transform( Belt_P2, Belt_P2, transformer( frame_type{C2 - G2, 0} ));
    geometry::transform( Belt_P1, Belt_P1, transformer( mark1, mark2 ));
    geometry::transform( Belt_P2, Belt_P2, transformer( mark2, mark1 ));

    value_type dist1 = 0, dist2 = 0;
    for (std::size_t i = 0; i != Belt_P1.size(); ++i)
        dist1 = std::max(dist1, distance(Belt_P1_copy[i] + C1, Belt_P1[i] + G2));
    for (std::size_t i = 0; i != Belt_P2.size(); ++i)
        dist2 = std::max(dist2, distance(Belt_P2_copy[i] + C2, Belt_P2[i] + G1));
    // version raccourcie + modifiée
    

    // %%%%%%%%%% Cas obj1 %%%%%%%%%%
    value_type calc;
    if (dist1 < 1e-12)
        calc = dt_defaut;
    else
        calc = ((d - lambda) / 2) * (dt_defaut / dist1);
    value_type dt12 = std::min(dt_defaut, calc);
    // %%%%%%%%%% Fin %%%%%%%%%%
    
    // %%%%%%%%%% Cas obj2 %%%%%%%%%%
    if (dist2 < 1e-12)
        calc = dt_defaut;
    else
        calc = ((d - lambda) / 2) * (dt_defaut / dist2);
    value_type dt21 = std::min(dt_defaut, calc);
    // %%%%%%%%%% Fin %%%%%%%%%%

    // assert( std::min(dt12, dt21) > 0 );
    if (std::min(dt12, dt21) < 0)
    {
        std::cout << d << " " << dist1 << " "  << " BUG DELTA T " << dist2;
        return 1e-4;
    }

    return std::min(dt12, dt21);
}

template <typename TDomain, typename TDetector>
typename TimeScaleManager<TDomain, TDetector>::value_type
TimeScaleManager<TDomain, TDetector>::delta_t_secu_fast(
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
        delta_t = - (( dist_secu - lambda ) / 2) / VRel;
        // delta_t = - ( dist_secu - lambda ) / VRel; // too much ?
    } else
    {
        // Collision impossible
        delta_t = std::numeric_limits<value_type>::max();
    }

    if (delta_t < 1e-12 ) std::cout
    << "BUG dtfast  " 
    << delta_t << " " 
    << dist_secu << " " 
    << lambda << " " 
    << VRel << " "
    << Vg1 << " "
    << Vg2
    <<std::endl; //DEBUG

    return delta_t;

}


}} // namespace floe::ope


#endif // OPE_TIME_SCALE_MANAGER_HPP