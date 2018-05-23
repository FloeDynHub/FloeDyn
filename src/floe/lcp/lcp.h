/*! \brief Definition of Linear Complementarity Problem
 *
 * \todo problem occuring when one want to compare integers of different signs: 'int' and 'std::size_t'. 
 *  perhaps modifying the dimension as an integers?
 *
 * \author Matthias Rabatel
 * \date April 2018
 * \copyright ...
 */

#ifndef DEF_FLOE_LCP
#define DEF_FLOE_LCP

#include <cassert>
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace floe { namespace lcp
{

/*! A LCP defined by M, q
 *
 * The problem is to find z that respect:
 *
 * w = Mz + q >= 0,
 *
 * z >= 0 and
 *
 * w . z = 0
 *
 */
template < typename T >
struct LCP
{
    //!< attributes:
    typedef boost::numeric::ublas::matrix<T> array_type;    //!< Type of array.
    typedef boost::numeric::ublas::vector<T> vector_type;   //!< Type of vector.

    std::size_t         dim;    //!< Dimension of the problem. (4 times the number of contacts)
    array_type          M;      //!< The matrix of the problem. Actually APS formulation \cite Moreau1988, Anitescu1997, Stewart2000.
    vector_type         q;      //!< The vector of the problem.
    vector_type         z;      //!< The solution of the LCP.
    int                 driving;//!< Index of the last found driving variable of the Lemke's algorithm.
                                //!< Useful whether Lemke's algorithm ends with a secondary ray.
    std::vector<int>    basis;  //!< The basis initialized with the canonical basis: 
                                //!< w variables are in {0,...,n-1} (basis variables)
                                //!< z variables are in {n,...,2n-1} (nonbasis variables).
                                //!< The solution is associated with a feasible basis \cite Cottle1992 p103.

    //!< functions:
    //! Constructor given the dimension of the problem.
    LCP( std::size_t n ) : dim(n), M(n, n), q(n), z(n,0.0), driving{-1}, basis(2*n+1) 
    {
        for (std::size_t i=0; i<2*n+1; ++i){
            basis[i] = i;
        }
    };

    //! Two constructors for the augmented LCP:
    LCP ( std::size_t n, array_type A );

    LCP ( std::size_t n, array_type A, vector_type d );

    /*! \fn T LCP_error( LCP<T> const& lcp, bool calc_w = true )
     * \brief Calculation of the LCP error
     *
     * Let E is lcp error, there exists two ways to compute E from the outputs, (z,w), of Lemke's algorithms:
     * E(z)     = sum( |z_i^-| + |( (Mz)_i + q_i )^-| + |z_i*(Mz+q)_i| )
     * E(z,w)   = sum( |z_i^-| +       |w_i^-|        +    |z_i*w_i|   )
     * 
     * In most cases, E(z) and E(z,w) are similar since w = Mz + q. However, one should be careful. First, in the case
     * where Lemke's algorithm terminates with a secondary ray, w = Mz + q + dz_0 with z_0 is still a basis. Second, it may
     * possible the iterative pivoting operations propagate numerical errors. After n iterations these errors could be lead to 
     * an unexpected w not anymore equal to Mz + q + dz_0, even when z_0 becomes a nonbasis variable. Thus, only E(z) is
     * selected for the LCP error computation.
     *
     * @param[in]   lcp     The linear complementarity to solve
     * @param[out]  err     A template parameter for the LCP error: E(z) or E(w,z).
     *
     * \todo One should remove the function LCP::LCP_error_detailed(LCP<T> const& lcp) after implementing the best solution
     * for LCP solvers.
     *
     */
    T LCP_error() const;
    
    /*! \fn boost::numeric::ublas::vector<T> LCP_error_detailed( LCP<T> const& lcp)
     * \brief Used for keeping track of the error localisation.
     */
    vector_type LCP_error_detailed() const;

    /*! \fn void pivoting(int row, int column)
     *  \brief Compute the pivoting operation on the matrix M with the pivot M(block, drive). 
     * 
     * @param[in] row       The row index corresponding to the blocking variable.
     * @param[in] column    The column index corresponding to the driving variable.
     */
    void pivoting(int row, int column);

    /*! \fn void multi_pivoting( matrix<T> subM, std::vector<int> idx_a)
     *  \brief Compute the pivoting operation by a sub-matrix non-singular
     *
     *  M square matrix:
     *  sub-Matrix = M_aa with a=idx_a the set of the z-basis.
     *  g=idx_g such as (idx_a,idx_g) is a partition of {0,...,dim-1}
     *  
     *        |   M_aa^-1    |      -M_aa^-1 M_ag        |
     *   M' = |--------------|---------------------------|
     *        | M_ga M_aa^-1 |  M_gg - M_ga M_aa^-1 M_ag |
     *  
     *        |     -M_aa^-1 q_a       |
     *   q' = |------------------------|
     *        | q_g - M_ga M_aa^-1 q_a |
     *  
     */
    void multi_pivoting(LCP<T> &lcp_orig, array_type subM, std::vector<int> idx_a);

    /*! \fn void reinit(LCP<T> &lcp_ori)
     *  \brief Used for lcp re-initialization 
     */
    void reinit(LCP<T> &lcp_ori);

    /*! \fn bool go_through_adj_cone()
     *  \brief Used for finding the solution modifying one and only one z-basis.
     *
     *  When the Lemke's algorithm ends with a secondary ray, one can go through an adjacent cone with principal pivoting operation:
     *  (pivot < Z0, Zr > with Z0 the variable of the additional covering vector and Zr the current unbounded variable) or 
     *  with pivoting operation of M with sub-matrix M_{\gamma \gamma} if it is nonsingular (\gamma is the set of the z-basis variables). 
     *  Advantages: These 2 pivoting techniques are similar if < Z0, Zr > is non-zero and M_{\alpha \alpha} is nonsingular (proof?). 
     *  Drawbacks: These techniques do not preserve the lexicographical ordering, so, cycling (i.e.: come down on a basis) is possible!!
     *  If both the pivot < Z0, Zr > is zero and M_{\alpha \alpha} is singular. 
     * 
     *  \warning If both the pivot < Z0, Zr > is zero and M_{\alpha \alpha} is singular, false is returned and nothing is done.
     *
     * \remark + in C++ 20 a std function for set difference is implementing. This could help for the pivoting operation with a sub-matrix.
     */
    bool go_through_adj_cone(LCP<T> &lcp_orig, const int Z0, const double tolerance);
}; // End Struct
  

}} // namespace floe::lcp

#endif // FLOE_LCP_LCP_HPP