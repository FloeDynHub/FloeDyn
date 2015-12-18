/*!
 * \file ope/generator_LCP_solver.hpp
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

namespace floe { namespace ope
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
    using value_type = T;

    GeneratorLCPSolver() : base_class() { base_class::epsilon = 0; }

private:


    //! Test LCP solution validity
    bool LCPtest(int compt, value_type EC, value_type born_EC, value_type Err, bool VRelNtest ) override;

};


template<typename T>
bool GeneratorLCPSolver<T>::LCPtest(int compt, value_type EC, value_type born_EC, value_type Err, bool VRelNtest ) {
    return !(EC > born_EC*(1+1e-2) || VRelNtest == 0);
}


}} // namespace floe::ope


#endif // OPE_GENERATOR_LCP_SOLVER_HPP