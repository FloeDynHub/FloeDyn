/*!
 * \brief Lemke algorithm with lexicographical ordering.
 * 
 * LCP solver based on an augmented LCP with a covering vector d = (1,...,1)^T.
 * Each pivoting iteration is performed with a single pivot < w_r, z_s >, where w_r is called
 * the blocking of the driving variable, denoted z_s.
 * The algorithm terminates with a secondary ray if the driving column contains only positive terms or with a solution
 * at hand if z_0 is the blocking variable.
 * 
 * \param[in] tolerance Represents the admissible tolerance for w negative. It is the same order as the
 * LCP error. 
 * 
 * \remark The tolerance is only used for increase the probability to take z_0 as pivot. Thus, in this case, after
 * the pivot operation, a solution is at hand satisfying: sum( |z_i^-| + |w_i^-| ) >= dim*tolerance, 
 * where dim is the dimension of the problem and the exponent v^- denote the negative coefficient of v.
 * If z_0 is not a part of possible pivots, then, the tolerance is not applied.
 *
 * More informations in \cite cottle1992 p265, 299, 352-357, Sect. 4.5, 4.9, 4.10.
 *
 * \author Matthias Rabatel
 * \date April 2018
 * \copyright ...
 */
#ifndef FLOE_LCP_SOLVER_LEXICOLEMKE_MR_H
#define FLOE_LCP_SOLVER_LEXICOLEMKE_MR_H

#include <vector>

#include "floe/lcp/lcp.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace floe { namespace lcp { namespace solver
{

using namespace boost::numeric::ublas;

template < typename T>
std::vector<int> lexicolemke_MR(double tolerance, LCP<T> &lcp, int itermax );

template<typename T>
std::vector<int> lcp_lexicolemke_MR( const double tolerance, const int itermax, const std::size_t dim, 
        matrix<T> &M, vector<T> &q, vector<T> &z, std::vector<int> &basis, int &driving  );

template<typename T>
void pivoting(matrix<T> &M, vector<T> &q, const int block, const int drive);

}}}

#endif