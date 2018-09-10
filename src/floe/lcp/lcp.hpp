/*! \brief Definition of Linear Complementarity Problem
 *
 * \author Matthias Rabatel
 * \date April 2018
 * \copyright ...
 */

#ifndef DEF_FLOE_LCP_HPP
#define DEF_FLOE_LCP_HPP

#include <assert.h>
#include <vector>
#include <algorithm>

#include "lcp.h"

// #include <boost/numeric/ublas/matrix.hpp>
// #include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>

using namespace boost::numeric::ublas;

namespace floe { namespace lcp
{

////////////////////////////////////////////////////////////
template < typename T >  
LCP<T>::LCP( std::size_t n, array_type A ) : dim(n), M(n, n+1), q(n), z(n,0.0), driving{-1}, basis(2*n+1) 
{
    for (std::size_t i=0; i<2*n+1; ++i){
        basis[i] = i;
    }

    assert(A.size1()==n);
    project(M,range(0,n),range(0,n)) = A;
    vector_type d(n,1);
    column(M,n) = d;
}

////////////////////////////////////////////////////////////
template < typename T >  
LCP<T>::LCP( std::size_t n, array_type A, vector_type d ) : dim(n), M(n, n+1), q(n), z(n,0.0), driving{-1}, basis(2*n+1)
{
    for (std::size_t i=0; i<2*n+1; ++i){
        basis[i] = i;
    }

    assert(A.size1()==n);
    project(M,range(0,n),range(0,n)) = A;
    column(M,n) = d;
}

////////////////////////////////////////////////////////////
template < typename T >
T LCP<T>::LCP_error() const
{
    // Calculating w = Mz + q
    auto w = prod(M, z) + q;

    // Error on w
    T w_err = 0;
    for ( T value : w ) if (value < 0) w_err -= value;

    // Error on z
    T z_err = 0;
    for ( T value : z ) if (value < 0) z_err -= value;

    // Error on z.w
    T zw_err = 0;
    auto itz = z.begin();
    auto itw = w.begin();
    for ( ; itz != z.end(); ++itz, ++itw )
    {
        if (*itz != 0 && *itw != 0) 
            zw_err += std::abs( (*itz) * (*itw) );
    }

    T err = w_err + z_err + zw_err;

    return err;
}

////////////////////////////////////////////////////////////
template < typename T >
vector<T> LCP<T>::LCP_error_detailed() const
{
    vector_type Vec_Err(3*dim,0);

    // Calculation: w = Az + q;
    vector_type w;
    w = prod(M, z) + q;

    // Calculation: zw = z^T w;
    vector_type zw(dim);
    for (std::size_t i=0; i<dim; ++i) {
        zw(i) = z(i) * w(i);
    }

    for (std::size_t i=0; i<dim; ++i) {
        Vec_Err(i) = zw(i); // energy part 
        if (z(i)<0) {Vec_Err(i+dim) = -z(i);} // impulse part
        if (w(i)<0) {Vec_Err(i+2*dim) = -w(i);} // relative velocities after contact
    }

    return Vec_Err;
}

////////////////////////////////////////////////////////////
template<typename T>
void LCP<T>::pivoting(const int block, const int drive)
{
    // WARNING: preferring a/b instead of 1/b*a
    // M is not necessarily a square matrix (augmented lcp case)
    // typedef typename LCP<T>::vector_type vector_type;
    assert(M(block,drive) !=0);
    
    std::size_t i, j;

    T pivot = M( block, drive );
    
    // M:
    vector_type d(dim,0.0);
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

////////////////////////////////////////////////////////////
template<typename T>
void LCP<T>::multi_pivoting( LCP<T> &lcp_orig, matrix<T> invSubM, std::vector<int> idx_a )
{
    std::vector<int> idx_g;
    std::vector<int>::iterator it;
    // /Warning: favourised the declaration of the indexes k, kc, kr inside the loop (more easy to debug!!)
    std::size_t k, kc, kr;

    for (k=0; k<dim; ++k) {
        it = std::find( idx_a.begin(), idx_a.end(), k );
        if (it==idx_a.end()) {
            idx_g.push_back(k);
        }
    }

    array_type   M_gg(idx_g.size(),idx_g.size()), M_ag(idx_a.size(),idx_g.size()), M_ga(idx_g.size(),idx_a.size()), 
        Mprime_gg(idx_g.size(),idx_g.size()), Mprime_ag(idx_a.size(),idx_g.size()), Mprime_ga(idx_g.size(),idx_a.size());

    vector_type  q_g(idx_g.size()), q_a(idx_a.size()), qprime_g(idx_g.size()), qprime_a(idx_a.size());

    // M_gg, q_g:
    for (kr=0; kr<idx_g.size(); ++kr) {
        q_g(kr) = lcp_orig.q(idx_g[kr]);
        for (kc=0; kc<idx_g.size(); ++kc) {
            M_gg(kr,kc) = lcp_orig.M(idx_g[kr],idx_g[kc]);
        }
    }

    // M_ag, q_a:
    for (kr=0; kr<idx_a.size(); ++kr) {
        q_a(kr) = lcp_orig.q(idx_a[kr]);
        for (kc=0; kc<idx_g.size(); ++kc) {
            M_ag(kr,kc) = lcp_orig.M(idx_a[kr],idx_g[kc]);
        }
    }

    // M_ga:
    for (kr=0; kr<idx_g.size(); ++kr) {
        for (kc=0; kc<idx_a.size(); ++kc) {
            M_ga(kr,kc) = lcp_orig.M(idx_g[kr],idx_a[kc]);
        }
    }

    // M'_gg, M'_ag, M'_ga:
    /*
     *  \Warning:   prod() no return a matrix or a vector, one need to transform to a matrix or a vector
     *              before use again prod().
     */
    Mprime_gg = M_gg - prod( M_ga , matrix<T>(prod(invSubM , M_ag)) );
    Mprime_ag = - prod( invSubM, M_ag );
    Mprime_ga = prod( M_ga, invSubM );
    qprime_a  = -prod( invSubM, q_a );
    qprime_g  = q_g - prod( M_ga , vector<T>(prod(invSubM, q_a)) ); 

    // reset to zero-matrix and zero-vector:
    M.resize(M.size1(),M.size2(),false);
    q.resize(q.size(),false);

    // M'_aa, q'_a:
    for (kr=0; kr<idx_a.size(); ++kr) {
        q(idx_a[kr]) = qprime_a(kr);
        for (kc=0; kc<idx_a.size(); ++kc) {
            M(idx_a[kr],idx_a[kc]) = invSubM(kr,kc);
        }
    }
    // M'_gg, q'_b:
    for (kr=0; kr<idx_g.size(); ++kr){
        q(idx_g[kr]) = qprime_g(kr);
        for (kc=0; kc<idx_g.size(); ++kc){
            M(idx_g[kr],idx_g[kc]) = Mprime_gg(kr,kc);
        }
    }
    // M_'ag
    for (kr=0; kr<idx_a.size(); ++kr) {
        for (kc=0; kc<idx_g.size(); ++kc) {
            M(idx_a[kr],idx_g[kc]) = Mprime_ag(kr,kc);
        }
    }
    // M_'ga
    for (kr=0; kr<idx_g.size(); ++kr) {
        for (kc=0; kc<idx_a.size(); ++kc) {
            M(idx_g[kr],idx_a[kc]) = Mprime_ga(kr,kc);
        }
    } 
    
    vector_type d(dim,1);
    column( M, dim ) = d;
}

////////////////////////////////////////////////////////////
template<typename T>
void LCP<T>::reinit(LCP<T> &lcp_ori)
{
    // works even if lcp_orig is an augmented LCP:
    project(M,range(0,dim),range(0,dim)) = project(lcp_ori.M, range(0,dim),range(0,dim));
    // if lcp is an augmented LCP:
    if (M.size1()+1==M.size2()) {
        vector_type d(dim,1); 
        column(M,dim)  = d;
    }
    q                                    = lcp_ori.q;
    z                                    = lcp_ori.z; 
    basis                                = lcp_ori.basis;
    driving                              = lcp_ori.driving;
}

////////////////////////////////////////////////////////////
template<typename T>
int LCP<T>::go_through_adj_cone( LCP<T> &lcp_orig, const int Z0, const double tolerance )
{
    int SR_status;
    std::vector<int>::iterator it;
    //! /warning favo
    std::size_t k, kc, kr;
    int block, drive;
    int n = dim; // definition of an integer equal to the dimension for comparison with an other integer.

    //preliminaries:
    it = std::find(basis.begin(), basis.begin()+dim-1, Z0);
    assert(*it==Z0);
    block = it-basis.begin();

    it = std::find(basis.begin()+dim, basis.end(), driving);
    assert(*it==driving);
    drive = it-basis.begin()-dim;
    
    if (M(block,drive) > tolerance/n) {// first use pivoting operation
        SR_status = 1;
        // std::cout << "inside go_through_adj_cone function with single pivoting!\n";
        // std::cout << "pivot: " << M(block,drive) << "\n\n";
        // std::cout << "block&drive: " << block << " | " << drive << "\n";
        this->pivoting( block, drive );
        // std::cout << "after pivot: \n";
        // after pivoting, one conserves only q and M without the driving
        // column (since that will be next d column):
        array_type M_temp(dim,dim+1);
        project( M_temp, range(0,dim), range(0,drive) ) = subrange( M, 0,dim , 0,drive );
        
        for (k=drive;k<dim;++k) {
            // std::cout << "column: " << k << "\n";
            column( M_temp, k) = column( M, k+1 );
            // updating basis and nonbasis
            basis[dim+k] = basis[dim+k+1];   
        }
        // std::cout << "driving: " << driving << "\n";
        basis[block] = driving;
        basis[2*dim] = Z0;

        vector_type d(dim,1);
        column( M_temp, dim ) = d;

        // std::cout << "basis: \n";
        // for (k=0;k<basis.size();++k) {
        //     std::cout << basis[k] << ", ";
        // } 
        // std::cout << "\n\n";

        M = M_temp;
    }    
    else {// second use direct inversing operation
        SR_status = 2;
        // std::cout << "inside go_through_adj_cone function with multi pivoting!\n";
        std::vector<int> idx_alpha_temp, idx_alpha, idx_beta, idx_rem;
        for (k=0;k<dim;++k) {
            if (z[k] > tolerance/n) {
                idx_alpha_temp.push_back(k);
            }
        }

        // equivalent to find z-basis.
        // Indeed if z_i = 0 is a z-basis then one can send back to
        // nonbasis.
        // Due to (APS) formulation, subM must not contain a
        // alpha-column alone (that means: (APS) formulation 
        // contains zero-diagonal entries corresponding to the alpha
        // coefficient of z the impact. If, after pivoting, there
        // exists a alpha-index in the basis, and there is no
        // beta-index associated with, then the submatrix will be
        // singular!!
        // Elimination of the alpha-index column.
        const int nb_pt = dim/4;
        assert(idx_alpha_temp.size()!=0);

        for (k=0;k<idx_alpha_temp.size();++k) {
            // removing "alpha" and "lambda" "-index": 
            if ( (idx_alpha_temp[k] > nb_pt-1) && (idx_alpha_temp[k] < 3*nb_pt) ) { 
                idx_beta.push_back( idx_alpha_temp[k]-nb_pt );
            }
        }

        for (k=0;k<idx_beta.size();++k) {
            idx_rem.push_back(idx_beta[k]/2 + 3*nb_pt);
        }

        for (k=0;k<idx_alpha_temp.size();++k) {
            if (idx_alpha_temp[k]<3*nb_pt) { // fullfilling with \lambda- and \beta-index:
                idx_alpha.push_back(idx_alpha_temp[k]);
            }
            else { // fullfilling with \alpha-index associated with beta-index:
                it = std::find(idx_rem.begin(),idx_rem.end(),idx_alpha_temp[k]);
                if (it!=idx_rem.end()) {
                    idx_alpha.push_back(idx_alpha_temp[k]);
                }
            }
        }
        
        assert(idx_alpha.size()!=0);
        matrix<T> subM(idx_alpha.size(),idx_alpha.size());
        
        for (kr=0;kr<idx_alpha.size();++kr){
            for (kc=0;kc<idx_alpha.size();++kc){
                subM(kr,kc) = lcp_orig.M(idx_alpha[kr],idx_alpha[kc]);
            }
        }

        // pivoting operation from a sub-matrix, only if the sub-matrix is nonsingular.
        // computation of the inverse of the sub-matrix: using the LU decomposition:
        // create a working copy of the subM:
        matrix<T> A(subM);
        // create a permutation matrix for the LU-factorization
        permutation_matrix<std::size_t> pm(A.size1());

        // perform LU-factorization
        int res = lu_factorize(A,pm);
        
        if( res != 0 ) {
            SR_status = 0;
            return SR_status; // sub-matrix is singular, the method is not feasible!
        }

        // create identity matrix of "inverse"
        matrix<T> invSubM(A.size1(),A.size2());
        invSubM.assign(identity_matrix<T>(A.size1()));

        // backsubstitute to get the inverse
        lu_substitute(A, pm, invSubM);

        // checking the quality of the inverse:
        // T N2 = norm_1( prod(subM, invSubM) - identity_matrix<T>(subM.size1()) );
        // assert(N2<1e-9);

        this->multi_pivoting( lcp_orig, invSubM, idx_alpha);

        // updating basis and nonbasis:
        for (k=0; k<2*dim+1; ++k){
            basis[k] = k;
        }
        for (k=0;k<idx_alpha.size();++k) {
            basis[idx_alpha[k]]     = idx_alpha[k]+dim;
            basis[idx_alpha[k]+dim] = idx_alpha[k];
        }
    }

    // this method is possible, that means either the pivot M(block,drive) is non zero or
    // the sub-matrix from the set of z-basis is non-singular.
    return SR_status;
}

}} // namespace floe::lcp

#endif // FLOE_LCP_LCP_HPP