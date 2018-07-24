/*!
 * \author Matthias Rabatel
 * \date May 2018
 * \copyright ...
 */
#ifndef FLOE_LCP_SOLVER_LEXICOLEMKE_MR_HPP
#define FLOE_LCP_SOLVER_LEXICOLEMKE_MR_HPP

#include <cassert>
#include <iostream>
#include <vector>

#include "floe/lcp/lcp.h"
#include "floe/lcp/solver/lexicolemke_MR.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace floe { namespace lcp { namespace solver
{

using namespace boost::numeric::ublas;

///////////////////////////////////////////////////////////////////////////////////////////////////
template < typename T>
std::vector<int> lexicolemke_MR(double tolerance, LCP<T>& lcp, int itermax)
{
    const double tol = tolerance;
    const int itmax = itermax;
    const std::size_t dim = lcp.dim;
    matrix<T> M = lcp.M;
    vector<T> q = lcp.q;
    vector<T> z = lcp.z;
    std::vector<int> bas = lcp.basis;
    int drive = lcp.driving;


    std::vector<int> error_status = lcp_lexicolemke_MR( tol, itmax, dim, M, q, z, bas, drive );

    lcp.M = M; lcp.q = q; lcp.z = z; lcp.basis = bas; lcp.driving = drive;

    return error_status;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
std::vector<int> lcp_lexicolemke_MR( const double tolerance, const int itermax, const std::size_t dim, 
        matrix<T> &M, vector<T> &q, vector<T> &z, std::vector<int> &basis, int &driving  )
{
    // itermax, by default: 1000
    int err{-1}, iter;
    int use_lexico_order = 0;
    std::vector<int> res(2,0);

    std::size_t i, j;
    int block, drive, entering, leaving, Z0_priority;
    double tol = tolerance/dim; // same order as the LCP error
    /* to prevent to remove all basis variable as candidate for pivoting (cause of strong numerical approximation)
     * one adjust the tolerance with the smallest value of q
     */
    std::vector<T> q_abs(3*dim/4,0);
    for (i=0;i<3*dim/4;++i) {q_abs[i] = std::abs(q[i]);}
    typename std::vector<T>::iterator tol_it = std::min_element( q_abs.begin(), q_abs.end() );
    T val_tol = *tol_it*1e-5;
    if (val_tol<tol) {tol = 1e-13;} // std::cout << "I am using the tol adjustement: " << *tol_it << " | " << tol << std::endl;}
    // End of tol adjustement

    std::vector<int>::iterator it, drive_it, Z0_priority_it;
    typename std::vector<T>::iterator min_it;    

    // trivial solution:
    int num_negative{0};
    for (i=0; i<dim; ++i) {
        if (q[i]>=0) {++num_negative;}
    }
    if ( num_negative == int(dim) ){
        z = q;
        res[0] = err; res[1] = use_lexico_order;
        return res;
    }   

    int Z0  = 2*dim; // artificial variable associated with the covering vector
    entering = Z0;
    
    // % w variables are in {0,...,dim-1}
    // % z variables are in {dim,...,2*dim-1}
    // % basis and nonbasis variables
    std::vector<int> bas, nonbas;
    for (i=0; i<dim; ++i) {
        bas.push_back(basis[i]); nonbas.push_back(basis[i+dim]);
    }
    nonbas.push_back(Z0);

    // Q = Id, matrix with lexicographically positive row built like a vector of vector for using 
    // lexicographical comparison:
    std::vector< std::vector<T> > Q(dim);
    for (i=0; i<dim; ++i) {
        std::vector<T> v(dim,0);
        v[i] = 1;
        Q[i] = v;
    }

    // Looking for zbar such as w = q + (d * zbar) >=0
    std::vector<T> x(dim,0);
    for (i=0; i<dim; ++i) {x[i] = q(i);}   
    min_it  = std::min_element( x.begin(), x.end() );
    block   = ( min_it - x.begin() );
    drive   = dim;
    
    leaving = bas[block];

    // keeping M and q before pivoting operations for checking of numerical error:
    decltype(M) M_orig = M;
    decltype(q) q_orig = q;

    // Pivoting operation on M and Q:
    T pivot = M( block, drive );
    // Q:
    // the blocking row:
    for (j=0; j<dim; ++j){
        Q[block][j] = -Q[block][j]/pivot; 
    }
    // the remaining of the matrix:
    for (i=0; i<dim; ++i) {
        for (j=0; j<dim; ++j) {             
            if (int(i)!=block) {
                Q[i][j] += M(i,drive)*Q[block][j];
            }
        }
    }

    // M:
    pivoting( M, q, block, drive );

    // updating of the nonbasis: 
    bas[block]      = entering;
    nonbas[drive]   = leaving;

    // variables used in the followinf loop:
    T val;                                  // minimul ratio test
    std::vector<T> d(dim,0);                // covering vector
    std::vector<int> E_block;               // index of negative values of d
    std::vector<T> r_t_tol, r_t;            // to compute the minimum ratio test
    std::vector<int> candidate_pivots_indx; // index of blocking variables satisfying the minimum
                                            // ratio test (MRT).
    std::vector<int> base_candidate;        // basis index for blocking satisfying the MRT
    std::size_t nb_candidate{0};            // length of candidate_pivots_indx
    std::vector<T> Q_tmp(dim,0.0);          // temporary vector for lexicographic comparison

    for (iter=1; iter<=itermax; ++iter) {

        if (leaving==Z0) {break;}
        else if ( leaving < int(dim) ) { // the blocking variable is a w
            entering = dim+leaving; // the complementarity variable of w is a z variable
        }
        else { // the blocking variable is a z
            entering = leaving - dim; // the complementarity variable of z is a w variable
        } 

        // updating of the covering vector:
        drive_it    = std::find( nonbas.begin(), nonbas.end(), entering);
        drive       = ( drive_it - nonbas.begin() );

        auto d = column( M, drive );

        // finding new blocking variable:
        for (i=0;i<dim;++i){
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
            r_t_tol[i]  = - ( q(E_block[i]) + tol ) / d(E_block[i]);
            r_t[i]      = -    q(E_block[i])        / d(E_block[i]);
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
                // // debug:
                // std::cout << "lexicographic ordering: \n";
                // std::cout << "Q: \n";
                // for (std::size_t iq=0;iq<nb_candidate; ++iq) {
                //     std::cout << "[ ";
                //     for (std::size_t jq=0;jq<dim;++jq) {
                //         std::cout << Q[candidate_pivots_indx[iq]][jq] << ", ";
                //     }
                //     std::cout << " ]\n";
                // }
                // std::cout << "\n\n";
                // // end of debug
                for (i=1;i<nb_candidate;++i){ // lexicographic minimum ratio
                    if (Q_tmp > Q[candidate_pivots_indx[i]]){
                        block = candidate_pivots_indx[i];
                        Q_tmp = Q[candidate_pivots_indx[i]];
                    }
                }
                // // debug:
                // std::cout << "Q_tmp: \n";
                // for (std::size_t iq=0;iq<dim; ++iq) {
                //     std::cout << Q_tmp[iq] << ", ";
                // }
                // std::cout << "\n\n";
                // std::cout << "dim | iter | block | drive: " << dim << " | " << iter << " | " 
                //     << block << " | " << drive << "\n";
                // std::cout << "previous pivot:" << pivot << "\n\n"; 
                // // end of debug
            }
            else { block = candidate_pivots_indx[0]; }
        }

        leaving = bas[block];

        // updating of [Mtilde,Qtilde]:
        assert(d[block]<0);
        // Pivoting operation on M and Q:
        T pivot = M( block, drive );
        // Q:
        // the blocking row:
        for (j=0; j<dim; ++j){
            Q[block][j] = -Q[block][j]/pivot; 
        }
        // the remaining of the matrix:
        for (i=0; i<dim; ++i) {
            for (j=0; j<dim; ++j) {             
                if (int(i)!=block) {
                    Q[i][j] += M(i,drive)*Q[block][j];
                }
            }
        }
        // M, q:
        pivoting( M, q, block, drive );

        // setting to 0 each qi such as i satisfies the min ratio test and
        // < wi,zi > is not the pivot! (since numerical error may lead to
        // negative qi, this could be dramatic in the following!)
        if ( ( Z0_priority_it == base_candidate.end() ) && (nb_candidate > 1) ) {
            for (i=0;i<candidate_pivots_indx.size();++i) {
                if (candidate_pivots_indx[i]!=block) {
                    q(candidate_pivots_indx[i]) = 0;
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
      
    if (iter >= itermax && leaving != Z0){err=1;}

    vector<T> w(dim,0.0);   
    // re-initialization of z:
    for (i=0; i<dim; ++i) {
        z[i] = 0.0;
    }

    basis[Z0] = nonbas[dim];
    for (i=0;i<dim;++i) {
        basis[i]                = bas[i];
        basis[i+dim]            = nonbas[i];

        if (bas[i]!=Z0) {
            if ( bas[i]<int(dim) ) {
                w[bas[i]]       = q[i];
            }
            else {
                z[bas[i]-dim]   = q[i];
            }    
        }  
    }

    // numerical error test:
    if (err==-1) { // Lemke's algorithme ends with a solution
        T N2 = norm_2( w - vector<T>( prod( subrange(M_orig,0,dim,0,dim) , z ) ) - q_orig );
        if (N2>tolerance) {err=0;}
    }

    // save the last driving variable in the case where Lemke's algorithm ends with a secondary ray:
    driving = entering;

    res[0] = err; res[1] = use_lexico_order;
    return res;
}

////////////////////////////////////////////////////////////
template<typename T>
void pivoting(matrix<T> &M, vector<T> &q, const int block, const int drive)
{
    // WARNING: preferring a/b instead of 1/b*a
    // M is not necessarily a square matrix (augmented lcp case)
    // typedef typename LCP<T>::vector_type vector_type;
    assert(M(block,drive) !=0);
    
    std::size_t i, j;
    std::size_t dim = M.size1();

    T pivot = M( block, drive );
    
    // M:
    vector<T> d(dim,0.0);
    for (i=0;i<dim;++i) {
        d(i) = M(i,drive);
    }

    // the driving column and the blocking row:
    column( M, drive ) = column( M, drive ) / pivot;
    row(    M, block ) = -row(   M, block ) / pivot;
    M( block, drive )  = 1.0 / pivot;

    // the remaining of the matrix:
    for (i=0; i<dim; ++i) {
        for (j=0; j<M.size2(); ++j) {           // in the case where LCP(q,M) is an augmented LCP
                                                // (M contains the covering vector d)
            if (int(i)!=block && int(j)!=drive) {
                M(i,j) += d(i)*M(block,j);
            }
        }
    }
    
    // q:
    q(block) = - q(block) / pivot;
    for (i=0; i<dim; ++i) {
        if (int(i)!=block) {
                q(i) += d(i)*q(block);
            }
    }
}

}}} // namespace floe::lcp::solver

#endif // FLOE_LCP_SOLVER_LEXICOLEMKE_MR_HPP
