/*!
 * \file ope/LCP_solver.hpp
 * \brief LCP solver
 * \author Quentin Jouet
 */

#ifndef OPE_LCP_SOLVER_HPP
#define OPE_LCP_SOLVER_HPP

#include "floe/ope/LCP_solver.h"

#include "floe/lcp/solver/lexicolemke.hpp"
#include "floe/lcp/solver/lemke_eigen.hpp"
#include "floe/lcp/builder/graph_to_lcp.hpp"
#include "floe/collision/contact_graph.hpp" // boost::edges (todo : should be hidden)
#include <algorithm>


namespace floe { namespace ope
{


template<typename T>
bool LCPSolver<T>::solve( lcp_type& lcp ) {
        using namespace floe::lcp::solver;
        return (lemke(lcp) || lexicolemke(lcp));
    }


template<typename T>
template<typename TContactGraph>
std::array<vector<typename LCPSolver<T>::value_type>, 2>
LCPSolver<T>::solve( TContactGraph& graph, bool& success ) {
    using namespace floe::lcp::solver;

    floe::lcp::builder::GraphLCP<value_type, decltype(graph)> graph_lcp( graph );
    auto lcp_orig = graph_lcp.getLCP();
    auto lcp = lcp_orig;

    std::vector<std::vector<int>> MatLCP{
        {1,1},{2,1},{3,1},
        {0,1},{1,1},{2,1},{3,1},
        {0,1},{1,1},{2,1},{3,1},
        {0,2},{1,2},{2,2},{3,2},
        {0,2},{1,2},{2,2},{3,2},
        {0,3},{1,3},{2,3},{3,3},
    };

    // %%%%%%%%%%%%%%%%%%%%%%%%
    // % phase de compression %
    // %%%%%%%%%%%%%%%%%%%%%%%%

    unsigned int comptchgt{0};
    bool solved{0};
    vector<value_type> Solc(3 * graph_lcp.nb_floes), floe_impulses(graph_lcp.nb_floes, 0);

    decltype(lcp.z) best_z;
    value_type best_Err = std::numeric_limits<value_type>::max();

    while (!solved)
    {
        if (comptchgt >= MatLCP.size()) // passer le contact
        {
            success = 0;
            return {{graph_lcp.W, floe_impulses}};
        }

        bool solver_success{0};
        switch (MatLCP[comptchgt][0])
        {
            case 0:
                lcp = random_perturbation(lcp, 1e-10);
                break;
            case 1:
                solver_success = lemke(lcp);
                break;
            case 2:
                solver_success = lexicolemke(lcp);
                break;
            case 3:
                // solver_success = iterlemke(lcp); // don't have Iterlemke...
                break;
                // matlab : zc = IterLemke(Aorigin, Qcorigin, 1e-11, best.zc);
        }
        comptchgt++;

        if (!solver_success)
            continue;

        if (std::any_of(lcp.z.begin(), lcp.z.end(), is_nan<double>))
            continue;

        lcp_orig.z = lcp.z;

         // As-t-on une meilleure solution ?
        auto Err = LCP_error(lcp_orig);
        if (Err < best_Err || best_Err == std::numeric_limits<value_type>::max())
        {
            best_z = lcp.z;
            best_Err = Err;
        }

         // On teste toujours la meilleure solution trouvÃ©e
        lcp_orig.z = best_z;
        Err = best_Err;
        
         // Solution correspondante
        Solc = calcSolc(graph_lcp, lcp_orig);

        if (std::any_of(Solc.begin(), Solc.end(), is_nan<double>))
            continue;

         // Energie cinetique, Erreur LCP & Vit rel Normale :
        auto ECc = calcEc(Solc, graph_lcp.M, graph_lcp.W);
        Err = LCP_error(lcp_orig);
        auto vitrelnormtest = VRelNtest(prod(trans(graph_lcp.J), Solc), graph);
        solved = LCPtest(MatLCP[comptchgt][1],ECc,1,Err,vitrelnormtest);
        // if (solved && Err > 1) std::cout << "LCPErr : " << Err << std::endl;
        // if (solved) std::cout << "*" << MatLCP[comptchgt][1];
    }

    // %%%%%%%%%%%%%%%%%%%%%%%%%%
    // % phase de decompression %
    // %%%%%%%%%%%%%%%%%%%%%%%%%%
    vector<value_type> Sold(3 * graph_lcp.nb_floes);
    if ( epsilon != 0 )
    {
        lcp_type lcp_d_orig = lcp = graph_lcp.getLCP_d(lcp_orig, Solc, epsilon);
        value_type born_sup = graph_lcp.born_sup_d(lcp_orig, epsilon);
        // if (born_sup > 1) std::cout << "*****BORN SUP***** " << born_sup << " ";

        comptchgt = 0;
        solved = 0;
        best_Err = std::numeric_limits<value_type>::max();

        while (!solved)
        {
            if (comptchgt >= MatLCP.size()) // passer le contact
            {
                success = 0;
                floe_impulses = graph_lcp.impulse_vector(lcp_orig, epsilon);
                return {{(1 + epsilon) * Solc - epsilon * graph_lcp.W, floe_impulses}};
            }

            bool solver_success{0};
            switch (MatLCP[comptchgt][0])
            {
                case 0:
                    lcp = random_perturbation(lcp, 1e-10);
                    break;
                    // [A,Qc] = divmat(A,Qc,Iheart);
                case 1:
                    solver_success = lemke(lcp);
                    break;
                case 2:
                    solver_success = lexicolemke(lcp);
                    break;
                case 3:
                    // solver_success = lexicolemke(lcp); // don't have Iterlemke yet...
                    break;
                    // zc = IterLemke(Aorigin, Qcorigin, 1e-11, best.zc);
            }
            comptchgt++;

            if (!solver_success)
                continue;

            if (std::any_of(lcp.z.begin(), lcp.z.end(), is_nan<double>))
                continue;

            lcp_d_orig.z = lcp.z;

             // As-t-on une meilleure solution ?
            auto Err = LCP_error(lcp_d_orig);
            if (Err < best_Err || best_Err == std::numeric_limits<value_type>::max())
            {
                best_z = lcp.z;
                best_Err = Err;
            }

             // On teste toujours la meilleure solution trouvée
            lcp_d_orig.z = best_z;
            Err = best_Err;
            
             // Solution correspondante
            Sold = calcSold(graph_lcp, lcp_orig, lcp_d_orig, Solc);

            if (std::any_of(Sold.begin(), Sold.end(), is_nan<double>))
                continue; // std::cout << "******NAN******";

             // Energie cinetique, Erreur LCP & Vit rel Normale :
            auto ECd = calcEc(Sold, graph_lcp.M, graph_lcp.W);
            Err = LCP_error(lcp_d_orig);
            auto vitrelnormtest = VRelNtest(prod(trans(graph_lcp.J), Sold), graph);
            solved = LCPtest(MatLCP[comptchgt][1], ECd, 1 + born_sup, Err, vitrelnormtest);

            //solution alternative pour la conservation de l'EC
            if (solved)
            {
                if (ECd > 1)
                {
                    // std::cout << "ECD>1 => " << ECd << " ";
                    Sold = (1 + epsilon) * Solc - epsilon * graph_lcp.W; // return this instead of Sold
                    floe_impulses = graph_lcp.impulse_vector(lcp_orig, epsilon);
                } else {
                    // Impulse calculation
                    floe_impulses = graph_lcp.impulse_vector(lcp_orig, lcp_d_orig, epsilon);
                }
            }
        }
    } else {
        success = 1;
        floe_impulses = graph_lcp.impulse_vector(lcp_orig);
        return {{Solc, floe_impulses}};
    }

    success = 1;
    return {{Sold, floe_impulses}};

}


template<typename T>
bool LCPSolver<T>::LCPtest(int compt, value_type EC, value_type born_EC, value_type Err, bool VRelNtest ){
    bool resp = 1;
        if (compt == 1)
        {
            if (EC > born_EC*(1+1e-4) || Err > 1e-11 || VRelNtest == 0)
                resp = 0;
        }
        else if (compt == 2)
        {
            if (EC > born_EC * (1+1e-4) || Err > 1e-8 || VRelNtest == 0)
                resp = 0;
        }
        else if (compt == 3)
        {
            if (EC > born_EC*(1+1e-2) || VRelNtest == 0)// || Err > 1)
                resp = 0;
        }
        return resp;
}


template<typename T>
template<typename Tmat, typename Tvect>
typename LCPSolver<T>::value_type 
LCPSolver<T>::calcEc(const Tvect& S, const Tmat& M, const Tvect& w){
    return inner_prod(prod(S, M), S) / inner_prod(prod(w, M), w);
}


template<typename T>
template<typename TGraphLCP>
vector<typename LCPSolver<T>::value_type>
LCPSolver<T>::calcSolc(TGraphLCP& graph_lcp, LCPSolver<T>::lcp_type& lcp)
{   
    const std::size_t m = graph_lcp.nb_contacts;
    return graph_lcp.W + prod(
        graph_lcp.invM,
        prod(graph_lcp.J, subrange(lcp.z, 0, m)) + prod(graph_lcp.D, subrange(lcp.z, m, 3*m))
    );
}

template<typename T>
template<typename TGraphLCP>
vector<typename LCPSolver<T>::value_type>
LCPSolver<T>::calcSold(TGraphLCP& graph_lcp, lcp_type& lcp_c, lcp_type& lcp_d, vector<value_type> Solc )
{   
    const std::size_t m = graph_lcp.nb_contacts;
    vector<value_type> ezc = epsilon * subrange(lcp_c.z, 0, m);
    return Solc + prod(
        graph_lcp.invM,
        prod(graph_lcp.J, subrange(lcp_d.z, 0, m) + ezc) + prod(graph_lcp.D, subrange(lcp_d.z, m, 3*m))
    );
}


template<typename T>
typename LCPSolver<T>::lcp_type
LCPSolver<T>::random_perturbation(lcp_type& lcp, value_type max){
    // version randomly modifying all non-zeros values
    for (auto it1 = lcp.A.begin1(); it1 != lcp.A.end1(); ++it1) {
        for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
            if (*it2 != 0)
                *it2 += random_real(max);
        }
    }
    /* // Matlab version (divmat), modifying only one value in lcp.A and one in lcp.q (worse results)
    auto m = lcp.dim / 4;
    std::size_t n1 = std::rand() % (3 * m), n2 = std::rand() % (3 * m); // Iheart variable (=3*m) in Matlab code
    lcp.A(n1,n2) *= (1 + random_real(max));
    lcp.q(n1) *= (1 + random_real(max));
    */
    return lcp; 
}


template<typename T>
template<typename TContactGraph>
bool LCPSolver<T>::VRelNtest(const vector<value_type>& V, const TContactGraph& graph){
    size_t contact_id = 0;
    for ( auto const& edge : make_iterator_range( boost::edges( graph ) ) )
    {
        // Foreach contact between this 2 floes ...
        for ( auto const& contact : graph[edge] )
        {
            if (V[contact_id] < 0)
            {
                value_type delta = V[contact_id] * 10; // 10 = DT_DEFAULT // get dt_defaut ?
                if (delta > contact.dist / 50)
                    return 0;
            }
            ++contact_id;
        }
    }
    return 1;
}


template<typename T>
T LCPSolver<T>::random_real(T max)
{
    return m_uniform_distribution(m_random_generator) * max;
}


}} // namespace floe::ope


#endif // OPE_LCP_SOLVER_HPP