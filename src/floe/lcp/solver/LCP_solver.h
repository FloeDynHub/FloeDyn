/*!
 * \file lcp/solver/LCP_solver2.h
 * \brief LCP solver
 * \author Quentin Jouet
 */

#ifndef OPE_LCP_SOLVER_H
#define OPE_LCP_SOLVER_H


#include "floe/lcp/lcp.hpp"
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/blas.hpp>
#include <random>
#include <tuple>

#include <iostream> // debug


namespace floe { namespace lcp { namespace solver
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
    using real_type = T;

    LCPSolver(real_type epsilon) : epsilon{epsilon}, m_random_generator{}, m_uniform_distribution{-1, 1},
        m_nb_solvers{1}, m_solver_stats(3, m_nb_solvers, 0) {}

    // ~LCPSolver(){
        // std::cout << "LCP solver stats : " << std::endl;
        // std::cout << m_solver_stats << std::endl;
        // for (auto i = 0; i< m_solver_stats.size1(); ++i )
        // {
        //     int tot = 0;
        //     for (auto j = 0; j< m_solver_stats.size2(); ++j ) tot += m_solver_stats(i,j);
        //     std::cout << tot << " ";
        // }
        // std::cout << std::endl;
    // }

    //! Solve LCP
    bool solve( lcp_type& lcp );
    template<typename TContactGraph>
    std::array<vector<real_type>, 2> solve( TContactGraph& graph, bool& success  );

    int nb_solver_run{0}; // test (nb call run_solver() in step)
    double chrono_solver{0.0}; // test perf
    double max_chrono_solver{0.0}; // test perf

protected:

    real_type epsilon; //!< energy restitution coeff
    std::default_random_engine m_random_generator;
    std::uniform_real_distribution<real_type> m_uniform_distribution;
    int m_nb_solvers;
    
    matrix<int> m_solver_stats; // test solver stats

    //! Random small perturbation of LCP
    lcp_type random_perturbation(lcp_type& lcp, real_type max);
    void random_perturbation2(lcp_type& lcp, real_type max);
    real_type random_real(real_type max);
    bool run_solver(lcp_type& lcp, int id);


    //! Test LCP solution validity
    virtual bool LCPtest(int compt, real_type EC, real_type born_EC, real_type Err, bool VRelNtest );

    //! Compute normalized Kinetic Energy
    template<typename Tmat, typename Tvect>
    real_type calcEc(const Tvect& S, const Tmat& M, const Tvect& w);

    //! Compute solution of compression phase
    template<typename TGraphLCP>
    vector<real_type>
    calcSolc(TGraphLCP& graph_lcp, lcp_type& lcp);

    //! Compute solution of decompression phase
    template<typename TGraphLCP>
    vector<real_type>
    calcSold(TGraphLCP& graph_lcp, lcp_type& lcp_c, lcp_type& lcp_d, vector<real_type> Solc);

    //! Normal relative speed test
    template<typename TContactGraph>
    bool VRelNtest(const vector<real_type>& V, const TContactGraph& graph);

};


template<typename T>
inline bool is_nan(const T t){
    return (t != t);
}


}}} // namespace floe::lcp::solver


#endif // OPE_LCP_SOLVER_H