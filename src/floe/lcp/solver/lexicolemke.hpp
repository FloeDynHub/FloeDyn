/* Siconos-Numerics, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/*!
 * \file floe/lcp/solver/lexicolemke.hpp
 * \brief LCP Solver using Lemke algorithm with lexicographical ordering.
 * \author Roland Denis
 */

#ifndef FLOE_LCP_SOLVER_LEXICOLEMKE_HPP
#define FLOE_LCP_SOLVER_LEXICOLEMKE_HPP

#include "floe/lcp/lcp.hpp"

namespace floe { namespace lcp { namespace solver
{

/*! Frontend to the lexicolemke LCP solver.
 * 
 * \tparam T    Fundamental type.
 * \param lcp   The linear complementarity problem to solve.
 * \return      true if the solver successed.
 */
template < typename T>
bool lexicolemke( floe::lcp::LCP<T>& lcp );

/*! LCP solver using Lemke algorithm with lexicographical ordering to avoid degenerate cases
 *
 * Original from SICONOS projet (see licence above).
 *
 * \param dim   Dimension of the LCP(M,q).
 * \param M     Pointer to M.
 * \param q     Pointer to q.
 * \param[out]  zlem    Pointer to z in which the solution will be stored.
 * \param[out]  wlen    Pointer to w in which Mz+q will be stored.
 * \param[out]  info    Equals 0 if solver successed.
 */
void lcp_lexicolemke(int dim, const double * M, const double * q, double *zlem , int *info);

}}} // namespace floe::lcp::solver

#endif // FLOE_LCP_SOLVER_LEXICOLEMKE_HPP

