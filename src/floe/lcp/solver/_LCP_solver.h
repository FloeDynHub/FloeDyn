/*!
 * \file lcp/solver/LCP_solver.h
 * \brief LCP solver
 * \author Quentin Jouet
 */

#ifndef OPE_LCP_SOLVER_H
#define OPE_LCP_SOLVER_H


#include "floe/lcp/lcp.hpp"
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/blas.hpp>
#include <random>

#include <iostream> // debug

 // Siconos solver benchmark
#include "floe/lcp/solver/siconos_bench.hpp"
// Siconos solver benchmark

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

    LCPSolver() : epsilon{0.4}, m_random_generator{}, m_uniform_distribution{-1, 1} {} // todo : should epsilon be runtime parameter ?

    //! Solve LCP
    bool solve( lcp_type& lcp );
    template<typename TContactGraph>
    std::array<vector<real_type>, 2> solve( TContactGraph& graph, bool& success  );

protected:

    real_type epsilon; //!< energy restitution coeff
    std::default_random_engine m_random_generator;
    std::uniform_real_distribution<real_type> m_uniform_distribution;

     // Siconos solver benchmark
    SiconosBench sic_bench;
    // Siconos solver benchmark

    //! Random small perturbation of LCP
    lcp_type random_perturbation(lcp_type& lcp, real_type max);
    real_type random_real(real_type max);

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


template<typename T> // todo : move elsewhere
inline bool is_nan(const T t){
    return (t != t);
}


}}} // namespace floe::lcp::solver


#endif // OPE_LCP_SOLVER_H