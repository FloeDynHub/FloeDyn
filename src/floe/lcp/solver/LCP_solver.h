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
#include <assert.h>
#include <Eigen/SVD> // saving matrix


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
        m_nb_solvers{2}, tolerance{1e-7,1e-4,1e-5,2.5e-6,2.5e-5} {}//, m_solver_stats(3, m_nb_solvers, 0) {}

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
    std::array<vector<real_type>, 2> solve( TContactGraph& graph, bool& success, int lcp_failed_stats[] );

    // int nb_solver_run{0}; // test (nb call run_solver() in step)
    double chrono_solver{0.0}; // test perf
    double max_chrono_solver{0.0}; // test perf

protected:

    real_type epsilon; //!< energy restitution coeff
    std::default_random_engine m_random_generator;
    std::uniform_real_distribution<real_type> m_uniform_distribution;
    int m_nb_solvers;
    
    double tolerance[5]; // tolerance: 1/ good one for the global error on the LCP 
    // 2/ minimal one for acceptance for the global error on the LCP
    // 3/ and 4/ good tolerances for, respectively, the relative normal velocities and the kinetic energy after contact
    // 5/ minimal one for acceptance for the kinetic energy after contact.

    // matrix<int> m_solver_stats; // test solver stats

    //! Random small perturbation of LCP: see #5 and matlab file to change the random perturbation routines. 
    // Keeping only the case 2/ and 3/ (=reduction_via_perturbation) to only use on IterLemke and Lexico 
    // with a coef setted to 1e-8.
    lcp_type random_perturbation(lcp_type& lcp, real_type max);
    void random_perturbation2(lcp_type& lcp, real_type max);

    /* \fn void reduction_via_perturbation(lcp_type& lcp, real_type max)
       \brief Create a perturbed LCP in the forms of the reductible LCP (see issue 17)

       The idea of perturbing a matrix M by small positive ε along the main diagonal prior to solving LCP(q,M) 
       was discussed in the LCP literature as a possible strategy to convert the problem into a possibly easier one 
       (see [Adler and Verma, 2011] and the discussion about regularization in [Cottle et al. 1992], 5.6). 
       The key for doing it successfully is the ability to easily convert a solution to the perturbed problem into a solution 
       to the original problem. Over the years it had been shown that for several classes, but by no means all, this strategy is workable.
       For M in C-matrix (co-positive), the proposition:
       Given z in SOL(q, M (ε)), it is possible to find (in polynomial time) either z in SOL(q, M ) or a certificate for 
       SOL(q, M ) = \emptyset
    */ 
    void reduction_via_perturbation(lcp_type& lcp, double alpha);


    real_type random_real(real_type max);
    void run_solver(lcp_type& lcp, int id);

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
    bool Rel_Norm_Vel_test(const vector<real_type>& V, const TContactGraph& graph);

    //! Saving M and q from dealt with LCP(M,q) (solved and unsolved) for further analysis
    //! Return a boolean to prevent the maximum capacity to store (ex: 50 000 LCP ~ 250 Mo)
    bool saving_LCP_in_hdf5(lcp_type& lcp, bool solved, int test_idx, int solver_idx, const vector<real_type>& Err,
    int w_fail, int min_how_is_solved );
    
    //! Saving information about the source of the error: 100 for LCP error, 20 for increase of Kinetic Energy,
    //! 3 for relative normale velocity that may cause an interpenetration. Ex: 123 correspond to all of sources.
    int which_failure(const vector<real_type>& Err, bool Is_pos_rel_norm_vel );

    int Is_solved(lcp_type& lcp, bool Is_pos_rel_norm_vel );

};


template<typename T>
inline bool is_nan(const T t){
    return (t != t);
}


}}} // namespace floe::lcp::solver


#endif // OPE_LCP_SOLVER_H