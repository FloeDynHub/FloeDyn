#ifndef FLOE_LCP_SOLVER_LEMKE_EIGEN_HPP
#define FLOE_LCP_SOLVER_LEMKE_EIGEN_HPP

#include <cmath>
#include <iostream>
#include <vector>

#include <stdexcept>
#undef eigen_assert
#define eigen_assert(x) \
  if (!(x)) { throw (std::runtime_error("Eigen assert fail")); }

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "unsupported/Eigen/src/IterativeSolvers/Scaling.h"

#include "floe/lcp/lcp.hpp"

#ifndef isnan
# define isnan(x) \
  (sizeof (x) == sizeof (long double) ? isnan_ld (x) \
  : sizeof (x) == sizeof (double) ? isnan_d (x) \
  : isnan_f(x))
static inline int isnan_f  (float       _x) { return _x != _x; }
static inline int isnan_d  (double      _x) { return _x != _x; }
static inline int isnan_ld (long double _x) { return _x != _x; }
#endif

#ifndef isinf
# define isinf(x) \
  (sizeof (x) == sizeof (long double) ? isinf_ld (x) \
  : sizeof (x) == sizeof (double) ? isinf_d (x) \
  : isinf_f(x))
static inline int isinf_f(float       _x)
{ return !isnan (_x) && isnan (_x - _x); }
static inline int isinf_d(double      _x)
{ return !isnan (_x) && isnan (_x - _x); }
static inline int isinf_ld(long double _x)
{ return !isnan (_x) && isnan (_x - _x); }
#endif

namespace floe { namespace lcp { namespace solver {

int lcp_lemke(const Eigen::MatrixXd& _M, const Eigen::VectorXd& _q, Eigen::VectorXd & _z);

bool validate(const Eigen::MatrixXd& _M, const Eigen::VectorXd& _z, const Eigen::VectorXd& _q);

/*! Frontend to the lemke LCP solver.
 * 
 * \tparam T    Fundamental type.
 * \param lcp   The linear complementarity problem to solve.
 * \return      true if the solver successed.
 */
template < typename T>
bool lemke( floe::lcp::LCP<T>& lcp );


}}} // namespace floe::lcp::solver

#endif // FLOE_LCP_SOLVER_LEMKE_EIGEN_HPP

