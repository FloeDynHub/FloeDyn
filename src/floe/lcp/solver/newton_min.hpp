/*!
 * \file floe/lcp/solver/newton_min.hpp
 * \brief Newton method based LCP solver.
 * \author Quentin Jouet
 */

#ifndef FLOE_LCP_SOLVER_NEWTON_MIN_HPP
#define FLOE_LCP_SOLVER_NEWTON_MIN_HPP

#include "floe/lcp/lcp.hpp"

namespace floe { namespace lcp { namespace solver
{

/*! Frontend to the Newton min LCP solver.
 * 
 * \tparam T    Fundamental type.
 * \param lcp   The linear complementarity problem to solve.
 * \return      true if the solver successed.
 */
template < typename T>
bool newton_min( floe::lcp::LCP<T>& lcp );

/*! Newton method based LCP solver
 * Original from SICONOS projet (see licence above).
 *
 * \param lcp   The linear complementarity problem to solve.
 * \param[out]  info    Equals 0 if solver successed.
 */
void lcp_newton_min(floe::lcp::LCP<double>& lcp, int *info);

}}} // namespace floe::lcp::solver

#endif // FLOE_LCP_SOLVER_NEWTON_MIN_HPP

