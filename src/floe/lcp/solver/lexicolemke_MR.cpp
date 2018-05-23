/*! \fn void lexicolemke_MR(std::list<int> basis)
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
 */
#ifndef FLOE_LCP_SOLVER_LEXICOLEMKE_MR_HPP
#define FLOE_LCP_SOLVER_LEXICOLEMKE_MR_HPP

#include <cassert>
#include <vector>

#include "floe/lcp/lcp.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace floe { namespace lcp { namespace solver
{

using namespace boost::numeric::ublas;

template < typename T>
void lexicolemke_MR(double tolerance, floe::lcp::LCP<T>& lcp, int itermax ){
// std::vector<int> lexicolemke_MR(double tolerance, floe::lcp::LCP<T>& lcp, int itermax ){

    // itermax, by default: 1000
    int err{-1}, iter;
    int use_lexico_order = 0;
    // std::vector<int> res(2,0);

    std::size_t i, j;
    int block, drive, entering, leaving, Z0_priority;
    double tol = tolerance/lcp.dim; // same order as the LCP error

    std::vector<int>::iterator it, drive_it, Z0_priority_it;
    typename std::vector<T>::iterator min_it;    

    // trivial solution:
    int num_negative{0};
    for (i=0; i<lcp.dim; ++i) {
        if (lcp.q[i]>=0) {++num_negative;}
    }
    if ( num_negative == int(lcp.dim) ){
        lcp.z = lcp.q;
        // return err;
        // res[0] = err; res[1] = use_lexico_order;
        // return res;
    }   

    int Z0  = 2*lcp.dim; // artificial variable associated with the covering vector
    entering = Z0;
    
    // % w variables are in {0,...,dim-1}
    // % z variables are in {dim,...,2*dim-1}
    // % basis and nonbasis variables
    std::vector<int> bas, nonbas;
    for (i=0; i<lcp.dim; ++i) {
        bas.push_back(lcp.basis[i]); nonbas.push_back(lcp.basis[i+lcp.dim]);
    }
    nonbas.push_back(Z0);

    // Q = Id, matrix with lexicographically positive row built like a vector of vector for using 
    // lexicographical comparison:
    std::vector< std::vector<T> > Q(lcp.dim);
    for (i=0; i<lcp.dim; ++i) {
        std::vector<T> v(lcp.dim,0);
        v[i] = 1;
        Q[i] = v;
    }
    
    // Looking for zbar such as w = q + (d * zbar) >=0
    std::vector<T> x(lcp.dim,0);
    for (i=0; i<lcp.dim; ++i) {x[i] = lcp.q(i);}   
    std::cout << "la\n";
    min_it  = std::min_element( x.begin(), x.end() );
    block   = ( min_it - x.begin() );
    drive   = lcp.dim;
    std::cout << "block and drive: " << block << " and " << drive << "\n";
    std::cout << "zbar: " << *min_it << "\n";
    
    leaving = bas[block];
    std::cout << "leaving: " << leaving << "\n";
    // keeping M and q before pivoting operations for checking of numerical error:
    decltype(lcp.M) M_orig = lcp.M;
    decltype(lcp.q) q_orig = lcp.q;

    // Pivoting operation on M and Q:
    T pivot = lcp.M( block, drive );
    // Q:
    // the blocking row:
    for (i=0; i<lcp.dim; ++i){
        Q[block][i] = -Q[block][i]/pivot; 
    }
    // the remaining of the matrix:
    for (i=0; i<lcp.dim; ++i) {
        for (j=0; j<lcp.dim+1; ++j) {             
            if (int(i)!=block) {
                Q[i][j] += lcp.M(i,drive)*Q[block][j];
            }
        }
    }
    // M:
    lcp.pivoting(block, drive);

    // updating of the nonbasis: 
    bas[block]      = entering;
    nonbas[drive]   = leaving;

    // variables used in the followinf loop:
    T val;                                  // minimul ratio test
    std::vector<T> d(lcp.dim,0);            // covering vector
    std::vector<int> E_block;               // index of negative values of d
    std::vector<T> r_t_tol, r_t;            // to compute the minimum ratio test
    std::vector<int> candidate_pivots_indx; // index of blocking variables satisfying the minimum
                                            // ratio test (MRT).
    std::vector<int> base_candidate;        // basis index for blocking satisfying the MRT
    std::size_t nb_candidate{0};            // length of candidate_pivots_indx
    std::vector<T> Q_tmp(lcp.dim,0.0);      // temporary vector for lexicographic comparison

    for (iter=1; iter<=itermax; ++iter) {

        if (leaving==Z0) {break;}
        else if ( leaving < int(lcp.dim) ) { // the blocking variable is a w
            entering = lcp.dim+leaving; // the complementarity variable of w is a z variable
        }
        else { // the blocking variable is a z
            entering = leaving - lcp.dim; // the complementarity variable of z is a w variable
        } 

        // updating of the covering vector:
        drive_it    = std::find( nonbas.begin(), nonbas.end(), entering);
        drive       = ( drive_it - nonbas.begin() );

        auto d = column( lcp.M, drive );

        // finding new blocking variable:
        for (i=0;i<lcp.dim;++i){
            if (d(i) < 0) {E_block.push_back(i);}
        }
        
        if (E_block.size()==0) {    // no new pivot, secondary ray! 
            err = 2;                // unbounded driving variable
            break;
        }

        // minimum ratio test:
        // tol represents the tolerance to allow w<0 with  w >= -tol. 
        // tol is of the same order as the LCP error (1e-7). Divided by the
        // dimension of the problem (n) to keep solution below the LCP error
        // (1e-7). Indeed: sum(z(z<0))+sum(W(W<0)) >= n*tol because the
        // solution contains at most n nonzero!
        r_t_tol.resize(E_block.size()); r_t.resize(E_block.size());
        for (i=0;i<E_block.size();++i) {
            r_t_tol[i]  = - ( lcp.q(E_block[i]) + tol ) / d(E_block[i]);
            r_t[i]      = -    lcp.q(E_block[i])        / d(E_block[i]);
        }

        min_it  = std::min_element( r_t_tol.begin(), r_t_tol.end() );
        val = *min_it;

        // taking pivot lower than val we ensure for all i, qi^(k+1) >= -tol.
        // Moreover, keeping the same, test r_t_tol we ensure that for all
        // step k, for all i, qi >= -1e-7 alors at the step k+1, for all i,
        // qi >= -1e-7
        for (i=0;i<E_block.size();++i) {
            if (r_t[i] <= val) {candidate_pivots_indx.push_back(E_block[i]);}
        }

        // First: test if Z0 is one a candidate for the blocking variable
        // It is better to try to catch Z0 as blocking variable with a
        // slightly tolerance than continue pivoting iterations without
        // approximation.
        base_candidate.resize(candidate_pivots_indx.size());
        for (i=0;i<candidate_pivots_indx.size();++i) {
            base_candidate[i] = bas[candidate_pivots_indx[i]];
        }
        Z0_priority_it = find( base_candidate.begin(), base_candidate.end(), Z0);

        if (Z0_priority_it != base_candidate.end()) {
            Z0_priority = std::distance(base_candidate.begin(), Z0_priority_it);
            block = candidate_pivots_indx[Z0_priority];
        }
        else { // Second: if Z0 is not a candidate, doing the pivoting operation without approximation
            min_it  = std::min_element( r_t.begin(), r_t.end() );
            val = *min_it;
            candidate_pivots_indx.clear(); // re-initializing of pivot candidates
            for (i=0;i<E_block.size();++i) {
                if (r_t[i] == val) {candidate_pivots_indx.push_back(E_block[i]);}
            }

            nb_candidate = candidate_pivots_indx.size();
            if (nb_candidate > 1) {// lexicographic comparison:
                ++use_lexico_order;
                Q_tmp = Q[candidate_pivots_indx[0]];
                block = candidate_pivots_indx[0];
                std::cout << "using LEXICOOO\n";
                for (i=1;i<nb_candidate;++i){
                    if (Q_tmp > Q[candidate_pivots_indx[i]]){
                        block = candidate_pivots_indx[i];
                        Q_tmp = Q[candidate_pivots_indx[i]];
                    }
                }
            }
            else { block = candidate_pivots_indx[0]; }
        }

        leaving = bas[block];

        std::cout << "iter: " << iter << entering << " and " << leaving << "\n";

        // updating of [Mtilde,Qtilde]:
        assert(d[block]<0);
        // Pivoting operation on M and Q:
        T pivot = lcp.M( block, drive );
        // Q:
        // the blocking row:
        for (i=0; i<lcp.dim; ++i){
            Q[block][j] = -Q[block][j]/pivot; 
        }
        // the remaining of the matrix:
        for (i=0; i<lcp.dim; ++i) {
            for (j=0; j<lcp.dim; ++j) {             
                if (int(i)!=block) {
                    Q[i][j] += lcp.M(i,drive)*Q[block][j];
                }
            }
        }
        // M, q:
        lcp.pivoting(block, drive);
        // setting to 0 each qi such as i satisfies the min ratio test and
        // < wi,zi > is not the pivot! (since numerical error may lead to
        // negative qi, this could be dramatic in the following!)
        if ( ( Z0_priority_it == base_candidate.end() ) && (nb_candidate > 1) ) {
            std::cout << "using approx toZ000\n";
            for (i=0;i<candidate_pivots_indx.size();++i) {
                if (candidate_pivots_indx[i]!=block) {
                    lcp.q(candidate_pivots_indx[i]) = 0;
                }
            }   
        }

        // updating of the nonbasis and the basis:
        bas[block]      = entering;
        nonbas[drive]   = leaving;

        // clearing vector:
        E_block.clear();
        candidate_pivots_indx.clear();

    }
      
    if (iter >= itermax && leaving != Z0){ err=1;}

    decltype(lcp.z) w(lcp.dim,0.0);   
    // re-initialization of z:
    for (i=0; i<lcp.dim; ++i) {
        lcp.z[i] = 0.0;
    }

    for (i=0;i<lcp.dim;++i) {
        lcp.basis[i]                    = bas[i];
        lcp.basis[i+lcp.dim]            = nonbas[i];

        if (bas[i]!=Z0) {
            if ( bas[i]<int(lcp.dim) ) {
                w[bas[i]]               = lcp.q[i];
            }
            else {
                lcp.z[bas[i]-lcp.dim]   = lcp.q[i];
            }    
        }  
    }

    std::cout << lcp.z << "\n";
    std::cout << "err: " << err << "\n";
    // numerical error test:
    if (err==-1) { // Lemke's algorithme ends with a solution
        T N2 = norm_2( w - prod(subrange(M_orig,0,lcp.dim,0,lcp.dim), lcp.z) - q_orig );
        std::cout << "pb??\n";
        if (N2>tolerance) { err=0; std::cout << "err numerical pb: " << err << "\n"; }
    }

    // save the last driving variable in the case where Lemke's algorithm ends with a secondary ray:
    lcp.driving = entering;
    std::cout << "pb there??\n";

    // return err;

    // res[0] = err; res[1] = use_lexico_order;
    // std::cout << "res: " << res[0] << res[1] << "\n";
    // return res;
}

}}} // namespace floe::lcp::solver

#endif // FLOE_LCP_SOLVER_LEXICOLEMKE_MR_HPP
