/*!
 * \file domain/time_scale_manager.hpp
 * \brief Time Scale Manager
 * \author Quentin Jouet
 */

#ifndef OPE_TIME_SCALE_MANAGER_HPP
#define OPE_TIME_SCALE_MANAGER_HPP

#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/multi_point.hpp"
#include "floe/geometry/frame/frame_transformers.hpp"

// #include "floe/collision/matlab/proximity_data.hpp"

#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace floe { namespace domain
{

/*! TimeScaleManager
 *
 * Operator for time step determination
 * i.e. finding maximal time step that ensures no floes interpenetration
 *
 */

template<typename TProxData>
class TimeScaleManager
{

public:
    using proximity_data_type = TProxData;
    using real_type = typename proximity_data_type::real_type;
    using floe_type = typename proximity_data_type::floe_type;
    using point_type = typename floe_type::point_type;
    using optim_type = typename proximity_data_type::optim_type;
    using frame_type = typename floe_type::frame_type;
    using multi_point_type = floe::geometry::MultiPoint<point_type>;
    using floe_interface_type = typename floe_type::floe_interface_type;
    using optim_interface_type = typename optim_type::optim_interface_type;

    TimeScaleManager() : m_prox_data{nullptr} {}

    /*!
     * Returns time step taking all floes into account, base on detector informations
     */
    template <typename TDomain>
    real_type delta_t_secu(TDomain* domain);

    inline void set_prox_data_ptr(proximity_data_type const* ptr) { m_prox_data = ptr; }

private:

    proximity_data_type const* m_prox_data;

    /*!
     * Returns maximal delta_t for 2 floes beeing close
     * Corresponds to gestion_temps() in Matlab code
     */
    real_type delta_t_secu(
        real_type dist_secu,
        real_type dist_opt,
        const floe_type& floe1,
        const floe_interface_type& floe2,
        const optim_type& optim1,
        const optim_interface_type& optim2,
        short I,
        real_type dt_default
    );

    /*!
     * Returns maximal delta_t for 2 floes not beeing close
     * Corresponds to gestion_temps_fast() in Matlab code
     */
    real_type delta_t_secu_fast(
        real_type dist_secu,
        const floe_type& floe1,
        const floe_interface_type& floe2,
        const optim_type& optim1,
        const optim_interface_type& optim2
    );

};

template <typename TDetector>
template <typename TDomain>
typename TimeScaleManager<TDetector>::real_type
TimeScaleManager<TDetector>::delta_t_secu(TDomain* domain)
{
    real_type dt_default = domain->default_time_step();
    
    if (m_prox_data->interpenetration())
    {
        // std::cout << "interpenetration detected_timescalemanagerL99" << std::endl;
        domain->set_time_step(domain->time_step() / 5);
        return domain->time_step();
    }
    // std::cout << "dt_default = " << dt_default << std::endl;
    real_type global_min_dt = dt_default;
    
    #pragma omp parallel for
    // std::cout << "size(i,j) = " << m_prox_data->size1() << ", " << m_prox_data->size2() << std::endl;
    for (std::size_t i = 0; i < m_prox_data->size1(); ++i)
    {
        for ( std::size_t j = i+ 1; j != m_prox_data->size2(); ++j ) //size1 or size2? 
        {
            real_type delta_t;
            if (m_prox_data->get_indic(i,j) == 0)
            {   
                delta_t = delta_t_secu_fast(
                    m_prox_data->get_dist_secu(i,j),
                    m_prox_data->get_floe(i), m_prox_data->get_floe_itf(j), m_prox_data->get_optim(i), m_prox_data->get_optim_itf(j));
            }
            else
            {
                delta_t = delta_t_secu(
                    m_prox_data->get_dist_secu(i,j), m_prox_data->get_dist_opt(i,j),
                    m_prox_data->get_floe(i), m_prox_data->get_floe_itf(j), m_prox_data->get_optim(i), m_prox_data->get_optim_itf(j),
                    m_prox_data->get_indic(i,j), dt_default);  
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
typename TimeScaleManager<TDetector>::real_type
TimeScaleManager<TDetector>::delta_t_secu(
    real_type dist_secu,
    real_type dist_opt,
    const floe_type& floe1,
    const floe_interface_type& floe2,
    const optim_type& optim1,
    const optim_interface_type& optim2,
    short I,
    real_type dt_default
){
    const point_type& Vg1 = floe1.state().speed, Vg2 = floe2.state().speed;
    const real_type& Vt1 = floe1.state().rot, Vt2 = floe2.state().rot;
    const point_type& C1 = optim1.global_disk().center, C2 = optim2.global_disk().center;
    const real_type& R1 = optim1.global_disk().radius, R2 = optim2.global_disk().radius;
    const point_type& G1 = floe1.state().pos, G2 = floe2.state().pos;
    const real_type& tau1 = optim1.tau(), tau2 = optim2.tau();
    const real_type& dc1 = optim1.cdist(), dc2 = optim2.cdist();
    // real_type lambda = std::min(dc1, dc2) / 20;

    // real_type d = std::max(dist_opt, dist_secu);
    // lambda = std::min(lambda, d / 20);

    // // Calcul du deplacement d un point par rapport aux reperes en t+dt.
    // // repere a l'instant t+dt_default :
    // frame_type mark1{G1 + dt_default * (Vg1 - Vg2), dt_default * Vt1};
    // frame_type mark2{G2 + dt_default * (Vg2 - Vg1), dt_default * Vt2};
    // // fin calcul repere

    // calcul ceinture de points: (pour la rotation)
    // auto D1 = (distance(C1, G1) + R1 - tau1);
    // auto D2 = (distance(C2,G2) + R2 - tau2);
    // // multi_point_type Belt_P1, Belt_P2;
    // for (int i=0; i<50; ++i)
    // {
    //     auto angle = 2 * i *  M_PI / 50;
    //     Belt_P1.push_back(point_type{D1 * cos(angle), D1 * sin(angle)});
    //     Belt_P2.push_back(point_type{D2 * cos(angle), D2 * sin(angle)});
    // }

    // multi_point_type Belt_P1_be, Belt_P2_be, Belt_P1_af, Belt_P2_af;

    // using namespace geometry::frame; // import transformer, itransformer
    // /*
    // // version matlab raccourcie
    // geometry::transform( Belt_P1, Belt_P1_af, transformer( frame_type{C1 - G1, 0} ));
    // geometry::transform( Belt_P2, Belt_P2_af, transformer( frame_type{C2 - G2, 0} ));
    // geometry::transform( Belt_P1_af, Belt_P1_af, transformer( mark1, mark2 ));
    // geometry::transform( Belt_P2_af, Belt_P2_af, transformer( mark2, mark1 ));

    // real_type dist1 = 0, dist2 = 0;
    // for (std::size_t i = 0; i != Belt_P1.size(); ++i)
    //     dist1 = std::max(dist1, distance(Belt_P1[i] + C1, Belt_P1_af[i] + G2));
    // for (std::size_t i = 0; i != Belt_P2.size(); ++i)
    //     dist2 = std::max(dist2, distance(Belt_P2[i] + C2, Belt_P2_af[i] + G1));
    // // version raccourcie
    // */

    // // Version 2
    // geometry::transform( Belt_P1, Belt_P1_be, transformer( frame_type{G1, 0} ));
    // geometry::transform( Belt_P2, Belt_P2_be, transformer( frame_type{G2, 0} ));
    // geometry::transform( Belt_P1, Belt_P1_af, transformer( mark1 ));
    // geometry::transform( Belt_P2, Belt_P2_af, transformer( mark2 ));
    
    // real_type dist1 = 0, dist2 = 0;
    // for (std::size_t i = 0; i != Belt_P1.size(); ++i)
    //     dist1 = std::max(dist1, distance(Belt_P1_be[i], Belt_P1_af[i]));
    // for (std::size_t i = 0; i != Belt_P2.size(); ++i)
    //     dist2 = std::max(dist2, distance(Belt_P2_be[i], Belt_P2_af[i]));
    // // END Version 2
    // 
    // // %%%%%%%%%% Cas obj1 %%%%%%%%%%
    // real_type calc;
    // if (dist1 < 1e-12)
    //     calc = dt_default;
    // else
    //     calc = ((d - lambda) / 2) * (dt_default / dist1);
    // real_type dt12 = std::min(dt_default, calc);
    // // %%%%%%%%%% Fin %%%%%%%%%%
    
    // // %%%%%%%%%% Cas obj2 %%%%%%%%%%
    // if (dist2 < 1e-12)
    //     calc = dt_default;
    // else
    //     calc = ((d - lambda) / 2) * (dt_default / dist2);
    // real_type dt21 = std::min(dt_default, calc);
    // // %%%%%%%%%% Fin %%%%%%%%%%

    // assert( std::min(dt12, dt21) > 0 );

    // return std::min(dt12, dt21);


    // [S.D.Brenner revamp]    
    real_type eta = std::min(dc1, dc2);
    point_type Axe = (G2 - G1) / distance(G1, G2);
    auto D1 = (distance(C1, G1) + R1 - tau1);
    auto D2 = (distance(C2,G2) + R2 - tau2);
    // Conservative guess for maximum point-point relative velocity (an abolutely max bound can also be written with ||Vg1|| + ||Vg2||, but I think this should be good enough)
    real_type VRel = abs(geometry::dot_product(Vg2-Vg1, Axe)) + abs(Vt1)*D1 + abs(Vt2)*D2; 

    real_type lambda = std::min(dc1, dc2) / 20;
    real_type d = std::max(dist_opt, dist_secu);
    real_type dijMax = std::max( abs(d-lambda), eta );
    
    real_type delta_t = dijMax/VRel; // in a contact, d_opt = eta (theoretically), so d>=eta

    // std::cout << "calc_slow; d/eta = " << d/eta << "; dt = " << delta_t << " s" << std::endl;
    // std::cout << "calc_slow; dt = " << delta_t << std::endl;
    // std::cout << "calc_slow; d = " << dist_secu << std::endl;
    return std::min( delta_t, dt_default);
    

}

template <typename TDetector>
typename TimeScaleManager<TDetector>::real_type
TimeScaleManager<TDetector>::delta_t_secu_fast(
        real_type dist_secu,
        const floe_type& floe1,
        const floe_interface_type& floe2,
        const optim_type& optim1,
        const optim_interface_type& optim2
){

    point_type Vg1 = floe1.state().speed, Vg2 = floe2.state().speed;
    point_type C1 = optim1.global_disk().center, C2 = optim2.global_disk().center;
    real_type R1 = optim1.global_disk().radius, R2 = optim2.global_disk().radius;
    real_type tau1 = optim1.tau(), tau2 = optim2.tau();
    point_type G1 = floe1.state().pos, G2 = floe2.state().pos;
    real_type dc1 = optim1.cdist(), dc2 = optim2.cdist();

    // Calcul de la marge lambda
    real_type lambda = std::min(dc1, dc2) / 20; // a.k.a., tol3 (see Rabatel thesis, pg. 88) 
    real_type eta = std::min(dc1, dc2);     // [S.D.Brenner]: dc1,dc2 are the values of "eta1, eta2" from the paper

    // Axe reliant chaque couple de floe
    // point_type Axe = (C2 - C1) / distance(C1, C2);
    point_type Axe = (G2 - G1) / distance(G1, G2);

    // Vitesse relative projetee sur l'axe
    real_type VRel = geometry::dot_product(Vg2 - Vg1, Axe);

    // Calcul du dt max
    real_type delta_t;
    real_type dij = distance(G1, G2)-(R1-tau1+distance(G1,C1))-(R2-tau2+distance(G2,C2)); //doing this here without changing global_disk might lead to issues... but ¯\_(ツ)_/¯

    // [ *I think* This doesn't work properly because the centre of the disks (C1,C2) is not necessarily the same as the centroid of the floe. Thus there's some difference in the effective translation velocity]

    if (VRel < 0)
    {
        // Collision possible
        // delta_t = - ( dist_secu - lambda ) / VRel;
        delta_t = - std::max(dij-lambda, (eta)/2) / VRel; // If dij<eta/2, use eta/2 to compute the timestep (floes might be actively colliding), otherwise use dij
    } else
    {
        // Collision impossible
        delta_t = std::numeric_limits<real_type>::max();
    }

    // std::cout << "calc_fast; d = " << dist_secu << " m; dt = " << delta_t << " s"<< std::endl;
    // std::cout << "calc_fast; dt = " << delta_t << std::endl;
    // std::cout << "X(:,n)="<< G1 <<"; Y(:,n)="<< G2 << "; n=n+1;" <<std::endl;
    return delta_t;

}


}} // namespace floe::domain


#endif // OPE_TIME_SCALE_MANAGER_HPP