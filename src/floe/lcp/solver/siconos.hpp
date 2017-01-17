/*!
 * \file floe/lcp/solver/siconos.hpp
 * \brief Interface for Siconos LCP Solvers.
 * \author Roland Denis
 */

#ifndef FLOE_LCP_SOLVER_SICONOS_HPP
#define FLOE_LCP_SOLVER_SICONOS_HPP

#include <siconos/LCP_Solvers.h>
#include <siconos/lcp_cst.h>

#include "floe/lcp/lcp.hpp"

namespace floe { namespace lcp { namespace solver
{

/*! Frontend to any solver from siconos.
 * 
 * \tparam T    Fundamental type.
 * \param lcp   The linear complementarity problem to solve.
 * \return      true if the solver successed.
 */
template < typename T>
bool siconos( floe::lcp::LCP<T>& lcp, LCP_SOLVER solver_id)
{
	auto A = createNumericsMatrixFromData(0, lcp.dim, lcp.dim, lcp.A.data().begin()); // first arg = storage (0-dense/1-sparse)
    auto LCP = new LinearComplementarityProblem();
    LCP->M = A;
    LCP->size = lcp.dim;
    LCP->q = lcp.q.data().begin();
    double* Zs{lcp.z.data().begin()};
    double* Ws{lcp.w.data().begin()};
    SolverOptions options;
    NumericsOptions global_options;
    global_options.verboseMode = 0;
    linearComplementarity_setDefaultSolverOptions(LCP, &options, solver_id);
    int info = linearComplementarity_driver(LCP, Zs,Ws, &options, &global_options);
    return (info == 0);
}


}}} // namespace floe::lcp::solver

#endif // FLOE_LCP_SOLVER_SICONOS_HPP

