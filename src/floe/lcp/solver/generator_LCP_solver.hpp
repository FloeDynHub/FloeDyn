/*!
 * \file lcp/solver/generator_LCP_solver.hpp
 * \brief LCP solver for generator
 * \author Quentin Jouet
 */

#ifndef OPE_GENERATOR_LCP_SOLVER_HPP
#define OPE_GENERATOR_LCP_SOLVER_HPP

#include "floe/lcp/solver/lexicolemke.hpp"
#include "floe/lcp/solver/lemke_eigen.hpp"
#include "floe/lcp/lcp.hpp"
#include "floe/lcp/builder/graph_to_lcp.hpp"
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/blas.hpp>
#include <algorithm>
#include <random>


#include <iostream> // debug

namespace floe { namespace lcp { namespace solver
{

/*! LCPSolver
 *
 * Operator for LCP solving
 * No decompression phase and minimal solution validation
 *
 */

using namespace boost::numeric::ublas;

template<typename T>
class GeneratorLCPSolver : public LCPSolver<T>
{

public:
    using base_class = LCPSolver<T>;
    using real_type = T;

    GeneratorLCPSolver(real_type epsilon) : base_class(epsilon) {}

// private:


    //! Test LCP solution validity
    // bool LCPtest(int compt, real_type EC, real_type born_EC, real_type Err, bool VRelNtest ) override;

};


// template<typename T>
// bool GeneratorLCPSolver<T>::LCPtest(int compt, real_type EC, real_type born_EC, real_type Err, bool VRelNtest ) {
//     return !(EC > born_EC*(1+1e-2) || VRelNtest == 0);


}}} // namespace floe::lcp::solver


#endif // OPE_GENERATOR_LCP_SOLVER_HPP