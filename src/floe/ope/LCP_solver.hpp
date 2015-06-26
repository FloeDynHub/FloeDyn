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


#include <iostream> // debug

namespace floe { namespace ope
{

/*! LCPSolver
 *
 * Operator for LCP solving
 *
 */

 using namespace boost::numeric::ublas;

class LCPSolver
{

public:
    using lcp_type = floe::lcp::LCP<double>;

    using value_type = double; // TODO get from trait

    const bool solve( lcp_type& lcp );
    template<typename TContactGraph>
    const vector<value_type> solve( TContactGraph& graph );

private:

    lcp_type random_perturbation(lcp_type& lcp, value_type max);

    // TODO put elsewhere ?
    value_type random_real(value_type max);

    // TODO declare all methods
    const bool LCPtest(int compt, value_type EC, value_type born_EC, value_type Err, bool VRelNtest );

    template<typename Tmat, typename Tvect>
    value_type calcEc(const Tvect& S, const Tmat& M, const Tvect& w);

    template<typename TGraphLCP>
    auto calcSol(TGraphLCP& graph_lcp, LCPSolver::lcp_type& lcp)
    -> decltype(graph_lcp.W);

    const bool VRelNtest(const vector<value_type>& V);

};


const bool LCPSolver::solve( lcp_type& lcp ) {
        using namespace floe::lcp::solver;
        // return lexicolemke(lcp);
        // return lexicolemke(random_perturbation(lcp, 1e-15));
        
        // const bool s1 = lemke(lcp);
        // if (s1) {std::cout << "error1 : " << LCP_error(lcp);}
        // lcp_type lcp2 = random_perturbation(lcp, 1e-15);
        // const bool s2 = lexicolemke(lcp2);
        // if (s2) {
        //     lcp.z = lcp2.z;
        //     std::cout << " error2 : " << LCP_error(lcp) << std::endl;
        // }
        // return s2;
        
        // return lemke(random_perturbation(lcp, 1e-20));
        return (lemke(lcp) || lexicolemke(lcp));
        // return (lexicolemke(random_perturbation(lcp, 1e-15)) ||
        //         lemke(random_perturbation(lcp, 1e-20)));
    }


template<typename TContactGraph>
const vector<typename LCPSolver::value_type>
LCPSolver::solve( TContactGraph& graph ) {
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
    int comptchgt{0};
    bool solved{0};
    const vector<value_type> Solc;

    decltype(lcp.z) best_z;
    value_type best_Err = std::numeric_limits<value_type>::max();

    while (!solved)
    {
        comptchgt++;
        if (comptchgt > MatLCP.size()) // passer le contact
            return graph_lcp.W; // TODO be sure this is the correct return

        bool success{0};
        switch (MatLCP[comptchgt][0])
        {
            case 0:
                lcp = random_perturbation(lcp, 1e-10); //TODO compare matlab "divmat"
                break;
                // [A,Qc] = divmat(A,Qc,Iheart);
            case 1:
                success = lemke(lcp);
                break;
            case 2:
                success = lexicolemke(lcp);
                break;
            case 3:
                success = lexicolemke(lcp); // don't have Iterlemke...
                break;
                // zc = IterLemke(Aorigin, Qcorigin, 1e-11, best.zc);
        }

        if (!success)
            continue;

        lcp_orig.z = lcp.z;

         // As-t-on une meilleure solution ?
        auto Err = LCP_error(lcp_orig);
        if (Err < best_Err)
        {
            best_z = lcp.z;
            best_Err = Err;
        }

         // On teste toujours la meilleure solution trouvÃ©e
        lcp_orig.z = best_z;
        Err = best_Err;
        
         // Solution correspondante
        auto Solc = calcSol(graph_lcp, lcp_orig);
        
         // Energie cinetique, Erreur LCP & Vit rel Normale :
        auto ECc = calcEc(Solc, graph_lcp.M, graph_lcp.W); // TODO be sure this is the correct 3rd arg
        Err = LCP_error(lcp_orig);
        // vitrelnormtest = VRelNtest(trans(graph_lcp.J) * Solc, dist, Mc); // TODO matlab version
        auto vitrelnormtest = VRelNtest(prod(trans(graph_lcp.J), Solc));
        solved = LCPtest(MatLCP[comptchgt][2],ECc,1,Err,vitrelnormtest);
    }
    return Solc;

}


const bool LCPSolver::LCPtest(int compt, value_type EC, value_type born_EC, value_type Err, bool VRelNtest ){
    bool resp = 1;
        if (compt == 1)
        {
            if (EC > born_EC*(1+1e-4) || abs(Err) > 1e-11 || VRelNtest == 0)
                resp = 0;
        }
        else if (compt == 2)
        {
            if (EC > born_EC * (1+1e-4) || abs(Err) > 1e-8 || VRelNtest == 0)
                resp = 0;
        }
        else if (compt == 3)
        {
            if (EC > born_EC*(1+1e-2) || VRelNtest == 0)
                resp = 0;
        }
        return resp;
}

template<typename Tmat, typename Tvect>
typename LCPSolver::value_type 
LCPSolver::calcEc(const Tvect& S, const Tmat& M, const Tvect& w){
    return inner_prod(prod(S, M), S) / inner_prod(prod(w, M), w); // TODO (0)? (0,0)?
}


template<typename TGraphLCP>
auto LCPSolver::calcSol(TGraphLCP& graph_lcp, LCPSolver::lcp_type& lcp)
-> decltype(graph_lcp.W)
{   
    std::size_t m = graph_lcp.J.size2();
    return graph_lcp.W + prod(
        graph_lcp.invM,
        prod(graph_lcp.J, subrange(lcp.z, 0, m)) + prod(graph_lcp.D, subrange(lcp.z, m, 3*m))
    );
}


typename LCPSolver::lcp_type
LCPSolver::random_perturbation(lcp_type& lcp, value_type max){
    for (auto i = 0; i < lcp.A.size1(); i++) {
        for (auto j = 0; j < lcp.A.size2(); j++) {
            lcp.A(i, j) = lcp.A(i, j) + random_real(max);
        }
    }
    return lcp; 
}


const bool LCPSolver::VRelNtest(const vector<value_type>& V){
    // TODO matlab version (more flexible)
    for (auto& v: V)
    {
        if (v < 0)
            return 0;
    }
    return 1;
}


// TODO put elsewhere ?
typename LCPSolver::value_type 
LCPSolver::random_real(value_type max)
{
    return (value_type)(std::rand() - RAND_MAX / 2) / RAND_MAX * max;
}


}} // namespace floe::ope


#endif // OPE_LCP_SOLVER_HPP