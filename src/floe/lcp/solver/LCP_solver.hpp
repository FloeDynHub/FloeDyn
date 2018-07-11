/*!
 * \file lcp/solver/LCP_solver.hpp
 * \brief LCP solver
 * \author Quentin Jouet
 */

#ifndef OPE_LCP_SOLVER_HPP
#define OPE_LCP_SOLVER_HPP
// #include <mpi.h>
#include "floe/lcp/solver/LCP_solver.h"

#include "floe/lcp/solver/lexicolemke.hpp"
#include "floe/lcp/solver/lexicolemke_MR.hpp"

#include "floe/lcp/solver/lemke_eigen.hpp"

#include "floe/lcp/builder/graph_to_lcp.hpp"
#include "floe/collision/contact_graph.hpp" // boost::edges (todo : should be hidden)
#include <algorithm>

#include <chrono>

using namespace boost::numeric::ublas;


namespace floe { namespace lcp { namespace solver
{

template<typename T>
template<typename TContactGraph>
std::array<vector<typename LCPSolver<T>::real_type>, 2>
LCPSolver<T>::solve( TContactGraph& graph, bool& success, int lcp_failed_stats[] ) {

    floe::lcp::builder::GraphLCP<real_type, decltype(graph)> graph_lcp( graph );
    auto lcp_orig = graph_lcp.getLCP();

    T                       best_err        = std::numeric_limits<T>::max();
    decltype(lcp_orig.z)    best_z;
    decltype(lcp_orig.M)    perturb_M       = lcp_orig.M;

    std::size_t i, j;

    T min_lcp = std::numeric_limits<T>::max();
    const std::size_t sd = 3*lcp_orig.dim/4; // size of the Delassus Matrix
    for (i=0; i<sd; ++i) {
        for (j=0; j<sd; ++j) {
            if (lcp_orig.M(i,j)!=0 && std::abs(lcp_orig.M(i,j))<min_lcp) {
                min_lcp = lcp_orig.M(i,j);
            }
        }
    }
    const T alpha = min_lcp * coef_perturb; // coefficient for the perturbation

    lcp_type lcp_a(lcp_orig.dim,lcp_orig.M);
    lcp_a.q = lcp_orig.q;

    // variables for storing LCP:
    // static bool is_full_storage  = false;
    // static bool bool_save_solved = true;
    // int w_fail{0};

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPRESSION PHASE: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vector<real_type> Solc(3 * graph_lcp.nb_floes), floe_impulses(graph_lcp.nb_floes, 0);

    // variables:
    bool solved=false;

    bool SR_status{0}, RP_status{0};
    int  last_status{0};             // -1 for RP_status, 0 for neither, 1 for SR_status.
    int  itermax=1000;
    std::vector<int> error_status(2,0);

    bool use_lexico_ordering{0};     // "true" if the lexicographic ordering has been used during the Lemke's algorithm.

    int count_attempt{0}, count_SR{0}, count_SR_failed{0}, count_RP{0};
    const int Z0  = 2*lcp_a.dim;    // artificial variable associated with the covering vector
    
    while (!solved && count_attempt<=ite_max_attempt) {

        error_status = lexicolemke_MR(tolerance, lcp_a, itermax);

        // Always comparing to the orignal one: (lcp.M could be perturbed)
        lcp_orig.z = lcp_a.z;

        T LCP_err = lcp_orig.LCP_error();

        if (LCP_err < best_err) {
            best_z = lcp_a.z;
            best_err = LCP_err;
        }

        // accurate solution is found
        if (LCP_err<=tolerance) {  
            // corresponding solution:
            Solc = calcSolc(graph_lcp, lcp_orig);

            last_status = SR_status-RP_status;
            solved = true; SR_status = 0; RP_status = 0;
        }
        // accurate solution not found
        else {    
            switch (error_status[0]) { // error_status[0] = -1: Lemke's algo ends with a solution
                case 0:
                    // std::cout << "numerical error propagation\n"; // go to matrix perturbation
                    RP_status = 1; SR_status = 0;
                    break;
                case 1:
                    // std::cout << "need for greater iteration number\n"; // go to matrix perturbation
                    if (itermax<8000) {itermax *= 2;}
                    RP_status = 1; SR_status = 0;
                    break;
                case 2:
                    // std::cout << "secondary ray, go through an adjacent cone\n"; // go through an adjacent cone
                    RP_status = 0; SR_status = 1;
                    if (count_SR > 2) {SR_status = 0; RP_status = 1;}
                    break;
            }
        }

        if (SR_status) { // secondary ray, go through an adjacent cone

            bool    is_done = lcp_a.go_through_adj_cone( lcp_orig, Z0, tolerance );

            if (is_done) {
                ++count_SR;
            }
            else { // the method consisting to go through an adjacent cone is not feasible!
                ++count_SR_failed;
                // std::cout << "SR failed, perturbation is required\n";
                RP_status = 1; SR_status = 0;
            }
        }

        if (RP_status) {
            // no solution has been found. One tries to perturb the LCP using addition of alpha Id:
            lcp_a.reinit(lcp_orig);

            reduction_via_perturbation( lcp_a.dim, perturb_M, alpha );
            project(lcp_a.M,range(0,lcp_a.dim),range(0,lcp_a.dim)) = perturb_M;
            ++count_RP;
        }

        use_lexico_ordering = (error_status[1]!=0)? 1:0; // if, at least, once the lexicographic ordering is used
                                                         // during the Lemke's algo, one store this information
                                                         // just to know.
        ++count_attempt;
    }

    // Mat
    // Saving data on LCP:
    // if (!is_full_storage){
    //     if ( ( !solved ) || ( bool_save_solved ) ){
    //         std::cout << "An more: " << bool_save_solved << " solved LCP stored!\n";

    //         lcp_orig.z                  = best_z;
    //         Solc                        = calcSolc(graph_lcp, lcp_orig);
    //         bool Is_pos_rel_norm_vel    = Rel_Norm_Vel_test(prod(trans(graph_lcp.J), Solc), graph);
    //         vector<real_type> err_det   = lcp_orig.LCP_error_detailed();
    //         w_fail                      = which_failure( err_det, Is_pos_rel_norm_vel );   
    //         last_status                 = SR_status-RP_status;   

    //         is_full_storage = saving_LCP_in_hdf5( lcp_orig, solved, count_attempt, count_RP, count_SR, count_SR_failed, 
    //             last_status, use_lexico_ordering, best_err, w_fail );
    //         bool_save_solved = false;
    //     }
    // }
    // End saving data on LCP
    // EndMat

    if (!solved) {
        // std::cout << "An unsolved LCP there!\n";
        // std::cout << "With a LCP error:" << best_err << "\n";
        // bool_save_solved     = true;
        lcp_failed_stats[0] += 1;

        success = 0;
        return {{graph_lcp.W, floe_impulses}};
    }

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DECOMPRESSION PHASE: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vector<real_type> Sold(3 * graph_lcp.nb_floes);

    if ( epsilon != 0 ) {

        lcp_type lcp_d_orig = graph_lcp.getLCP_d(lcp_orig, Solc, epsilon);
        lcp_a.reinit(lcp_d_orig);
        perturb_M = lcp_d_orig.M;

        // variables:
        solved=false;
        itermax=1000;
        count_attempt=0;

        while (!solved && count_attempt<=ite_max_attempt) {

            error_status = lexicolemke_MR(tolerance, lcp_a, itermax);

            // Always comparing to the original one: (lcp.M could be perturbed)
            lcp_d_orig.z = lcp_a.z;
            T LCP_err = lcp_d_orig.LCP_error();

            // accurate solution is found
            if (LCP_err<=tolerance) {    
                Sold = calcSold(graph_lcp, lcp_orig, lcp_d_orig, Solc);
                auto ECd = calcEc(Sold, graph_lcp.M, graph_lcp.W);
                if (ECd > 1) {
                    // std::cout << "Oops I'm with an exceed of kinetic energy!\n";
                    lcp_failed_stats[2] += 1;
                    Sold = (1 + epsilon) * Solc - epsilon * graph_lcp.W; // return this instead of Sold
                    floe_impulses = graph_lcp.impulse_vector(lcp_orig, epsilon);
                } else {
                    // Impulse calculation
                    floe_impulses = graph_lcp.impulse_vector(lcp_orig, lcp_d_orig, epsilon);
                } 
                solved = true; SR_status = 0; RP_status = 0;
            }
            // accurate solution not found
            else {                          
                switch (error_status[0]) {
                    case 0:
                        // std::cout << "numerical error propagation\n"; // go to matrix perturbation
                        RP_status = 1; SR_status = 0;
                        break;
                    case 1:
                        // std::cout << "need for greater iteration number\n"; // go to matrix perturbation
                        if (itermax<8000) {itermax *= 2;}
                        RP_status = 1; SR_status = 0;
                        break;
                    case 2:
                        // std::cout << "secondary ray, go through an adjacent cone\n"; // go through an adjacent cone
                        SR_status = 0; RP_status = 1;
                        if (count_SR >= 2) {SR_status = 0; RP_status = 1;}
                        break;
                }
            }

            if (SR_status) { // secondary ray, go through an adjacent cone
                bool is_done = lcp_a.go_through_adj_cone( lcp_d_orig, Z0, tolerance );

                if (!is_done) { // the method consisting to go through an adjacent cone is not feasible!
                    RP_status = 1;
                }            
            }

            if (RP_status) {
                // no solution has been found. One tries to perturb the LCP using addition of alpha Id:
                lcp_a.reinit(lcp_d_orig);

                reduction_via_perturbation( lcp_a.dim, perturb_M, alpha );
                project(lcp_a.M,range(0,lcp_a.dim),range(0,lcp_a.dim)) = perturb_M;
            }

            ++count_attempt;
        }

        if (!solved) {
            lcp_failed_stats[1] += 1;

            success = 0;
            floe_impulses = graph_lcp.impulse_vector(lcp_orig, epsilon);
            /*! 
             *  \attention when, in the decompression phase, the LCP remains unsolved, we return the solution
             *             given by the linear combination of the velovies before and after the compression phase.
             */
            return {{(1 + epsilon) * Solc - epsilon * graph_lcp.W, floe_impulses}};
        } 
    }
    else {
        success = 1;
        floe_impulses = graph_lcp.impulse_vector(lcp_orig);
        return {{Solc, floe_impulses}};
    }

    success = 1;
    return {{Sold, floe_impulses}};
}

template<typename T>
template<typename Tmat, typename Tvect>
typename LCPSolver<T>::real_type 
LCPSolver<T>::calcEc(const Tvect& S, const Tmat& M, const Tvect& w)
{
    return inner_prod(prod(S, M), S) / inner_prod(prod(w, M), w);
}


template<typename T>
template<typename TGraphLCP>
vector<typename LCPSolver<T>::real_type>
LCPSolver<T>::calcSolc(TGraphLCP& graph_lcp, LCPSolver<T>::lcp_type& lcp)
{   
    const std::size_t m = graph_lcp.nb_contacts;
    return graph_lcp.W + prod(
        graph_lcp.invM,
        prod(graph_lcp.J, subrange(lcp.z, 0, m)) + prod(graph_lcp.D, subrange(lcp.z, m, 3*m))
    );
}

template<typename T>
template<typename TGraphLCP>
vector<typename LCPSolver<T>::real_type>
LCPSolver<T>::calcSold(TGraphLCP& graph_lcp, lcp_type& lcp_c, lcp_type& lcp_d, vector<real_type> Solc )
{   
    const std::size_t m = graph_lcp.nb_contacts;
    vector<real_type> ezc = epsilon * subrange(lcp_c.z, 0, m);
    return Solc + prod(
        graph_lcp.invM,
        prod(graph_lcp.J, subrange(lcp_d.z, 0, m) + ezc) + prod(graph_lcp.D, subrange(lcp_d.z, m, 3*m))
    );
}

template<typename T>
template<typename TContactGraph>
bool LCPSolver<T>::Rel_Norm_Vel_test(const vector<real_type>& V, const TContactGraph& graph)
{
    size_t contact_id = 0;
    for ( auto const& edge : make_iterator_range( boost::edges( graph ) ) )
    {
        // Foreach contact between this 2 floes ...
        for ( auto const& contact : graph[edge] )
        {
            if (V[contact_id] < 0)
            {
                real_type delta = - V[contact_id] * 10; // 10 = DT_DEFAULT // get dt_defaut ?
                if (delta > contact.dist / 50)
                    return false;
            }
            ++contact_id;
        }
    }
    return true;
}
  
template<typename T>
int LCPSolver<T>::which_failure( vector<real_type> Err, bool Is_pos_rel_norm_vel )
{
    std::size_t dim = Err.size()/3; // total dimension (3*4*nbc)
    std::size_t nbc = dim/4; // number of contact

    real_type EC_err(0), impulse_err(0);
    // global LCP error:
    for (std::size_t i=dim; i<2*dim; ++i) {
        if (Err(i)<0) {impulse_err += std::abs(Err(i));}
    }
    // error on the kinetic energy
    for (std::size_t i=0; i<3*nbc; ++i) {
        if (Err(i)>0) {EC_err += Err(i);}
    }

    int source=0;

    if (impulse_err>tolerance){source += 100;}
    if (EC_err>tolerance) {source += 20;}
    if (!Is_pos_rel_norm_vel){source += 3;}

    return source;
}

}}} // namespace floe::lcp::solver


#endif // OPE_LCP_SOLVER_HPP