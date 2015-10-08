/*!
 * \file ope/LCP_solver.hpp
 * \brief LCP solver
 * \author Quentin Jouet
 */

#ifndef OPE_LCP_SOLVER_HPP
#define OPE_LCP_SOLVER_HPP

#include "floe/lcp/solver/lexicolemke.hpp"
#include "floe/lcp/solver/lemke_eigen.hpp"
#include "floe/lcp/lcp.hpp"
#include "floe/lcp/builder/graph_to_lcp.hpp"
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/blas.hpp>
#include <algorithm>


#include <iostream> // debug

namespace floe { namespace ope
{

/*! LCPSolver
 *
 * Operator for LCP solving
 *
 */

using namespace boost::numeric::ublas;

template<typename T>
class LCPSolver
{

public:
    using lcp_type = floe::lcp::LCP<T>;
    using value_type = T;

    LCPSolver() : epsilon{0.4} {} // should epsilon be runtime parameter ?

    bool solve( lcp_type& lcp );
    template<typename TContactGraph>
    vector<value_type> solve( TContactGraph& graph, bool& success  );

private:

    value_type epsilon; // energy restitution coeff

    lcp_type random_perturbation(lcp_type& lcp, value_type max);

    // TODO put elsewhere
    value_type random_real(value_type max);

    bool LCPtest(int compt, value_type EC, value_type born_EC, value_type Err, bool VRelNtest );

    template<typename Tmat, typename Tvect>
    value_type calcEc(const Tvect& S, const Tmat& M, const Tvect& w);

    template<typename TGraphLCP>
    vector<value_type>
    calcSolc(TGraphLCP& graph_lcp, lcp_type& lcp);

    template<typename TGraphLCP>
    vector<value_type>
    calcSold(TGraphLCP& graph_lcp, lcp_type& lcp_c, lcp_type& lcp_d, vector<value_type> Solc);

    template<typename TContactGraph>
    bool VRelNtest(const vector<value_type>& V, const TContactGraph& graph);

};


template<typename T>
inline bool is_nan(const T t){
    return (t != t);
}


template<typename T>
bool LCPSolver<T>::solve( lcp_type& lcp ) {
        using namespace floe::lcp::solver;
        return (lemke(lcp) || lexicolemke(lcp));
    }


template<typename T>
template<typename TContactGraph>
vector<typename LCPSolver<T>::value_type>
LCPSolver<T>::solve( TContactGraph& graph, bool& success ) {
    using namespace floe::lcp::solver;

    floe::lcp::builder::GraphLCP<value_type, decltype(graph)> graph_lcp( graph );
    auto lcp = graph_lcp.getLCP();
    auto lcp_orig = lcp;

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
    vector<value_type> Solc(graph_lcp.J.size1());

    decltype(lcp.z) best_z;
    value_type best_Err = std::numeric_limits<value_type>::max();

    while (!solved)
    {
        if (comptchgt >= MatLCP.size()) // passer le contact
        {
            success = 0;
            return graph_lcp.W;
        }

        bool success{0};
        switch (MatLCP[comptchgt][0])
        {
            case 0:
                lcp = random_perturbation(lcp, 1e-10);
                break;
                // [A,Qc] = divmat(A,Qc,Iheart);
            case 1:
                success = lemke(lcp);
                break;
            case 2:
                success = lexicolemke(lcp);
                break;
            case 3:
                // success = iterlemke(lcp); // don't have Iterlemke...
                break;
                // zc = IterLemke(Aorigin, Qcorigin, 1e-11, best.zc);
        }
        comptchgt++;

        if (!success)
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
    }

    // %%%%%%%%%%%%%%%%%%%%%%%%%%
    // % phase de decompression %
    // %%%%%%%%%%%%%%%%%%%%%%%%%%
    vector<value_type> Sold(graph_lcp.J.size1());
    if ( epsilon != 0 )
    {
        lcp_type lcp_d_orig = lcp = graph_lcp.getLCP_d(lcp_orig, Solc, epsilon);
        value_type born_sup = graph_lcp.born_sup_d(lcp_orig, epsilon);

        comptchgt = 0;
        solved = 0;
        best_Err = std::numeric_limits<value_type>::max();

        while (!solved)
        {
            if (comptchgt >= MatLCP.size()) // passer le contact
            {
                success = 0;
                return (1 + epsilon) * Solc - epsilon * (graph_lcp.W);
            }

            bool success{0};
            switch (MatLCP[comptchgt][0])
            {
                case 0:
                    lcp = random_perturbation(lcp, 1e-10);
                    break;
                    // [A,Qc] = divmat(A,Qc,Iheart);
                case 1:
                    success = lemke(lcp);
                    break;
                case 2:
                    success = lexicolemke(lcp);
                    break;
                case 3:
                    // success = lexicolemke(lcp); // don't have Iterlemke yet...
                    break;
                    // zc = IterLemke(Aorigin, Qcorigin, 1e-11, best.zc);
            }
            comptchgt++;

            if (!success)
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

             // On teste toujours la meilleure solution trouvÃ©e
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
        }
    } else {
        return Solc;
    }

    success = 1;
    return Sold;

}


template<typename T>
bool LCPSolver<T>::LCPtest(int compt, value_type EC, value_type born_EC, value_type Err, bool VRelNtest ){
    bool resp = 1;
        if (compt == 1)
        {
            if (EC > born_EC*(1+1e-4) || std::abs(Err) > 1e-11 || VRelNtest == 0)
                resp = 0;
        }
        else if (compt == 2)
        {
            if (EC > born_EC * (1+1e-4) || std::abs(Err) > 1e-8 || VRelNtest == 0)
                resp = 0;
        }
        else if (compt == 3)
        {
            if (EC > born_EC*(1+1e-2) || VRelNtest == 0)
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
    const std::size_t m = graph_lcp.J.size2();
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
    const std::size_t m = graph_lcp.J.size2();
    vector<value_type> ezc = epsilon * subrange(lcp_c.z, 0, m);
    return Solc + prod(
        graph_lcp.invM,
        prod(graph_lcp.J, subrange(lcp_d.z, 0, m)) + prod(graph_lcp.D, subrange(lcp_d.z, m, 3*m)) + prod(graph_lcp.J, ezc)
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
    /* // Matlab version (divmat), modifying only one value in lcp.A and one in lcp.q (bad results)
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
                value_type delta = V[contact_id] * DT_DEFAULT; // get dt_defaut otherwise ?
                if (delta > contact.dist / 50)
                    return 0;
            }
            ++contact_id;
        }
    }
    return 1;
}


template<typename T>
typename LCPSolver<T>::value_type 
LCPSolver<T>::random_real(value_type max) // TODO put elsewhere
{
    return (value_type)(std::rand() - RAND_MAX / 2) / RAND_MAX * max;
}


}} // namespace floe::ope


#endif // OPE_LCP_SOLVER_HPP