/*!
 * \file lcp/solver/LCP_solver.h
 * \brief LCP solver
 * \author Quentin Jouet
 */

#ifndef OPE_LCP_SOLVER_H
#define OPE_LCP_SOLVER_H

#include "floe/lcp/lcp.h"
// #include <boost/numeric/ublas/vector_proxy.hpp>
// #include <boost/numeric/ublas/blas.hpp>

#include <iostream> // debug
#include <assert.h>
#include <Eigen/SVD> // saving lcp

// saving matrix when lcp solver failed for further analysing
#include "H5Cpp.h"                    

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

    /** A constructor
     *
     *  the tolerance is equal to 1e-7 for simulations without obstacles and 1e-6 for simulations with obstacles (due to 
     *  lighter weight and greater velocity).
     * 
     */
    LCPSolver(real_type epsilon) : epsilon{epsilon}, m_tolerance{1e-7}, m_coef_perturb{1e-8}, m_ite_max_attempt{10} {}

    template<typename TContactGraph>
    std::array<vector<real_type>, 2> solve( TContactGraph& graph, bool& success, int lcp_failed_stats[] );

    double chrono_solver{0.0}; // test perf
    double max_chrono_solver{0.0}; // test perf

protected:
    typedef boost::numeric::ublas::matrix<real_type> array_type;
    typedef boost::numeric::ublas::matrix<real_type> vector_type;

    real_type   epsilon;          /*!< energy restitution coefficient*/
    real_type   m_tolerance;        // tolerance consented for the LCP error
    real_type   m_coef_perturb;     // multiplier coefficient for the size of the perturbations (1e-9 seems to be an optimal value). 
                                // One could be taking from 1e-7 to 1e-10.
    int         m_ite_max_attempt;  // integer for the number of perturbations (5 seems to be a good compromise between 
                                // do not loose too much time and good succes rate). To increase the success rate, 
                                // one could increase up to 10.

    //! Compute normalized Kinetic Energy
    template<typename Tmat, typename Tvect>
    real_type calcEc(const Tvect& S, const Tmat& M, const Tvect& w);

    //! Compute solution of compression phase
    template<typename TGraphLCP>
    vector<real_type> calcSolc(TGraphLCP& graph_lcp, lcp_type& lcp);

    //! Compute solution of decompression phase
    template<typename TGraphLCP>
    vector<real_type> calcSold(TGraphLCP& graph_lcp, lcp_type& lcp_c, lcp_type& lcp_d, vector<real_type> Solc);

    //! Normal relative speed test
    template<typename TContactGraph>
    bool Rel_Norm_Vel_test(const vector<real_type>& V, const TContactGraph& graph);
    
    //! Saving information about the source of the error: 100 for LCP error, 20 for increase of Kinetic Energy,
    //! 3 for relative normale velocity that may cause an interpenetration. Ex: 123 correspond to all of sources.
    int which_failure( vector<real_type> Err, bool Is_pos_rel_norm_vel );

};

/*  \fn void reduction_via_perturbation(lcp_type& lcp, real_type max)
 *  \brief Create a perturbed LCP in the forms of the reductible LCP (see issue 17)
 *
 *  The idea of perturbing a matrix M by small positive ε along the main diagonal prior to solving LCP(q,M) 
 *  was discussed in the LCP literature as a possible strategy to convert the problem into a possibly easier one 
 *  (see \cite Adler2011 and the discussion about regularization in \cite cottle1992, 5.6). 
 *  The key for doing it successfully is the ability to easily convert a solution to the perturbed problem into a solution 
 *  to the original problem. Over the years it had been shown that for several classes, but by no means all, this strategy is workable.
 *  For M in C-matrix (co-positive), the proposition:
 *  Given z in SOL(q, M (ε)), it is possible to find (in polynomial time) either z in SOL(q, M ) or a certificate for 
 *  SOL(q, M ) = \emptyset
 */ 
template<typename T>
void reduction_via_perturbation(std::size_t dim , matrix<T> &M, T alpha);

/*  \fn saving_LCP_in_hdf5(lcp_type lcp, bool solved, int count_attempt, int count_RP, 
 *      int count_SR, int count_SR_failed, int last_status, int perturb_used, bool use_lexico_ordering, 
 *      real_type lcp_err, int w_fail)
 *  
 *  \brief  Saving M and q from dealt with LCP(M,q) (solved and unsolved) for further analysis
 *          Return a boolean to prevent the maximum capacity to store (ex: 50 000 LCP ~ 250 Mo).
 *
 *  \remark The storage used hdf5 file formulation and consists in the following kind of table:
 *          |     1     |     2      |     3      |   4   |          5         |   6    |      8      |        9       |
 *          |:---------:|:----------:|:----------:|:-----:|:------------------:|:------:|:-----------:|:--------------:|
 *          | LCP error | nb attempt | nb perturb | nb SR | nb adj cone failed | lexico | idx failure | last technique |
 *
 *          with: \e nb for number, \e SR for secondary ray, \e adj \e cone \e failed for the failure of the 
 *          method consisting in going through an adjacent cone, \e lexico is true if the lexicographic 
 *          ordering is used during at least one Lemke's algorithm, \e idx \e perturb for the index of the
 *          matrix perturbation (see LCPSolver::matrix_perturbation(const std::size_t dim , matrix &M, const double alpha, const int Idx_perturb)),
 *          \e idx \e failure for the source of the LCP error (see LCPSolver::which_failure( vector<real_type> Err, bool Is_pos_rel_norm_vel ))
 *          and \e last \e technique for the last SR or perturbation used before solving or last attempt done. 
 */
template<typename T>
bool saving_LCP_in_hdf5(floe::lcp::LCP<T> lcp, int m_ite_max_attempt, std::vector<double> stats_vec_lcp, bool solved, int w_fail);
// bool saving_LCP_in_hdf5(floe::lcp::LCP<T> lcp, int m_ite_max_attempt, bool solved, int count_attempt, int count_RP, 
//     int count_SR, int count_SR_failed, int last_status, bool use_lexico_ordering, T lcp_err, int w_fail);

}}} // namespace floe::lcp::solver


#endif // OPE_LCP_SOLVER_H