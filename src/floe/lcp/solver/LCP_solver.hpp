/*!
 * \file lcp/solver/LCP_solver.hpp
 * \brief LCP solver
 * \author Quentin Jouet
 */

#ifndef OPE_LCP_SOLVER_HPP
#define OPE_LCP_SOLVER_HPP
// #include <mpi.h>
#include "floe/lcp/solver/LCP_solver.h"

#include "floe/lcp/solver/lexicolemke_MR.hpp"

#include "floe/lcp/builder/graph_to_lcp.hpp"
#include "floe/collision/contact_graph.hpp" // boost::edges (todo : should be hidden)
#include <algorithm>

#include <chrono>

using namespace boost::numeric::ublas;


namespace floe { namespace lcp { namespace solver
{

template<typename T>
template<typename TContactGraph>
vector<typename LCPSolver<T>::real_type>
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
    const T alpha = min_lcp * m_coef_perturb; // coefficient for the perturbation

    lcp_type lcp_a(lcp_orig.dim,lcp_orig.M);
    lcp_a.q = lcp_orig.q;

    // variables for storing LCP:
    #ifdef LCPSTATS
        static bool is_full_storage     = false;
        static bool bool_save_solved    = true;
        int         w_fail              = 0;
    #endif

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPRESSION PHASE: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vector<real_type> Solc(3 * graph_lcp.nb_floes);

    // variables:
    bool solved=false;
    bool SR_status{0}, RP_status{0};
    int  itermax=1000;
    std::vector<int>    error_status(2,0);
    std::vector<double> stats_vec_lcp(2*m_ite_max_attempt,0);

    int count_attempt{1}, count_SR{0};
    const int Z0  = 2*lcp_a.dim;    // artificial variable associated with the covering vector
    
    while (!solved && count_attempt<=m_ite_max_attempt) {

        error_status = lexicolemke_MR(m_tolerance, lcp_a, itermax);

        // Always comparing to the orignal one: (lcp.M could be perturbed)
        lcp_orig.z = lcp_a.z;

        T LCP_err = lcp_orig.LCP_error();

        stats_vec_lcp[2*(count_attempt-1)]  = LCP_err;
        stats_vec_lcp[2*(count_attempt-1)+1]    = double(error_status[0])*10+double(error_status[1]); 

        if (LCP_err < best_err) {
            best_z = lcp_a.z;
            best_err = LCP_err;
        }

        // accurate solution is found
        if (LCP_err<=m_tolerance) {  
            // corresponding solution:
            Solc = calcSolc(graph_lcp, lcp_orig);

            solved = true; SR_status = 0; RP_status = 0;
        }
        // accurate solution not found
        else {    
            switch (error_status[0]) { // error_status[0] = -1: Lemke's algo ends with a solution
                case -1: // bad solution
                    RP_status = 1; SR_status = 0;
                    break;
                case -2: // bad solution
                    RP_status = 1; SR_status = 0;
                    break;
                case 3:
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
                    if (count_SR > 5) {
                        SR_status = 0; RP_status = 1; 
                        stats_vec_lcp[2*(count_attempt-1)+1]    = 10*stats_vec_lcp[2*(count_attempt-1)+1] + 3;
                    }
                    break;
            }
        }

        if (SR_status) { // secondary ray, go through an adjacent cone

            int    is_done = lcp_a.go_through_adj_cone( lcp_orig, Z0, m_tolerance ); // is_done = 0 for failed SR
            // is_done = 1 for successful SR with unique pivoting operation
            // is_done = 2 for successful SR with direct inversing operation

            if (is_done==0) { // the method consisting to go through an adjacent cone is not feasible!
                // std::cout << "SR failed, perturbation is required\n";
                RP_status = 1; SR_status = 0;
            } else {++count_SR;}
            stats_vec_lcp[2*(count_attempt-1)+1]    = 10*stats_vec_lcp[2*(count_attempt-1)+1] + double(is_done);
        }

        if (RP_status) {
            // no solution has been found. One tries to perturb the LCP using addition of alpha Id:
            lcp_a.reinit(lcp_orig);

            reduction_via_perturbation( lcp_a.dim, perturb_M, alpha );
            project(lcp_a.M,range(0,lcp_a.dim),range(0,lcp_a.dim)) = perturb_M;

            stats_vec_lcp[2*(count_attempt-1)+1]    = 10*stats_vec_lcp[2*(count_attempt-1)+1] + double(RP_status);
        }

        ++count_attempt;
    }

    // Mat
    // Saving data on LCP:
    #ifdef LCPSTATS
        if (!is_full_storage){
            if ( ( !solved ) || ( bool_save_solved ) ){
                // std::cout << "A more: " << bool_save_solved << " solved LCP stored!\n";

                lcp_orig.z                  = best_z;
                Solc                        = calcSolc(graph_lcp, lcp_orig);
                bool Is_pos_rel_norm_vel    = Rel_Norm_Vel_test(prod(trans(graph_lcp.J), Solc), graph);
                vector<real_type> err_det   = lcp_orig.LCP_error_detailed();
                w_fail                      = which_failure( err_det, Is_pos_rel_norm_vel );     

                is_full_storage = saving_LCP_in_hdf5( lcp_orig, m_ite_max_attempt, stats_vec_lcp, solved, w_fail );
                bool_save_solved = false;
            }
        }
    #endif
    // End saving data on LCP
    // EndMat

    if (!solved) {
        // std::cout << "An unsolved LCP there!\n";
        // std::cout << "With a LCP error:" << best_err << "\n";
        #ifdef LCPSTATS
            bool_save_solved     = true;
        #endif
            
        lcp_failed_stats[0] += 1;

        success = 0;
        return graph_lcp.W;
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
        count_attempt=1;
        count_SR=0;

        while (!solved && count_attempt<=m_ite_max_attempt) {

            error_status = lexicolemke_MR(m_tolerance, lcp_a, itermax);

            // Always comparing to the original one: (lcp.M could be perturbed)
            lcp_d_orig.z = lcp_a.z;
            T LCP_err = lcp_d_orig.LCP_error();

            // accurate solution is found
            if (LCP_err<=m_tolerance) {    
                Sold = calcSold(graph_lcp, lcp_orig, lcp_d_orig, Solc);
                auto ECd = calcEc(Sold, graph_lcp.M, graph_lcp.W);
                if (ECd > 1) {
                    // std::cout << "Oops I'm with an exceed of kinetic energy!\n";
                    lcp_failed_stats[2] += 1;
                    Sold = (1 + epsilon) * Solc - epsilon * graph_lcp.W; // return this instead of Sold
                    graph_lcp.compute_impulses(lcp_orig, epsilon);
                } else {
                    // Impulse calculation
                    graph_lcp.compute_impulses(lcp_orig, lcp_d_orig, epsilon);
                } 
                solved = true; SR_status = 0; RP_status = 0;
            }
            // accurate solution not found
            else {                          
                switch (error_status[0]) {
                    case -1:
                        RP_status = 1; SR_status = 0;
                        break;
                    case -2:
                        RP_status = 1; SR_status = 0;
                        break;    
                    case 3:
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
                        if (count_SR > 5) {SR_status = 0; RP_status = 1;}
                        break;
                }
            }

            if (SR_status) { // secondary ray, go through an adjacent cone
                int    is_done = lcp_a.go_through_adj_cone( lcp_orig, Z0, m_tolerance ); // is_done = 0 for failed SR
                // is_done = 1 for successful SR with unique pivoting operation
                // is_done = 2 for successful SR with direct inversing operation

                if (is_done==0) { // the method consisting to go through an adjacent cone is not feasible!
                    // std::cout << "SR failed, perturbation is required\n";
                    RP_status = 1; SR_status = 0;
                } else {++count_SR;}            
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
            graph_lcp.compute_impulses(lcp_orig, epsilon);
            /*! 
             *  \attention when, in the decompression phase, the LCP remains unsolved, we return the solution
             *             given by the linear combination of the velovies before and after the compression phase.
             */
            return (1 + epsilon) * Solc - epsilon * graph_lcp.W;
        } 
    }
    else {
        success = 1;
        graph_lcp.compute_impulses(lcp_orig);
        return Solc;
    }

    success = 1;
    return Sold;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
template<typename Tmat, typename Tvect>
typename LCPSolver<T>::real_type 
LCPSolver<T>::calcEc(const Tvect& S, const Tmat& M, const Tvect& w)
{
    return inner_prod(prod(S, M), S) / inner_prod(prod(w, M), w);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////
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
  
///////////////////////////////////////////////////////////////////////////////////////////////////////////  
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

    if (impulse_err>m_tolerance){source += 100;}
    if (EC_err>m_tolerance) {source += 20;}
    if (!Is_pos_rel_norm_vel){source += 3;}

    return source;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
void reduction_via_perturbation(std::size_t dim , matrix<T> &M, T alpha){

    const std::size_t size_Delassus = 3*dim/4;

    for (std::size_t i=0; i<size_Delassus; ++i) {
        for (std::size_t j = 0; j<size_Delassus; ++j) {
            if (i == j) {
                M(i,j) += alpha;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
bool LCPSolver<T>::saving_LCP_in_hdf5(floe::lcp::LCP<T> lcp, int m_ite_max_attempt, std::vector<double> stats_vec_lcp, bool solved, int w_fail)
{
    using namespace H5;                                 // proper way inside a function (to prevent extension to entire code)

    const H5std_string FILE_NAME( "io/outputs/LCP_stats.h5" );
    const H5std_string GROUP_NAME_I( "solved" ); // root group
    const H5std_string GROUP_NAME_II( "unsolved" ); // root group
    const H5std_string GROUP_NAME1( "M" ); 
    const H5std_string GROUP_NAME2( "q" );
    const H5std_string GROUP_NAME3( "z" );
    const H5std_string LCP_error( "LCP error" );
    const H5std_string Last_Memb( "Last LCP" ); // to prevent similar LCP
    const H5std_string Contact_Graph_Info( "Contact Graph Info" ); // to store information on the contact graph and the "while loop"
    // const H5std_string Idx_solver( "Which solver" ); // Information on which solver, 
    // how many random perturbations are used before to compute solution, the index in the h5 file and
    // the source of the LCP error (see which_failure)

    const hsize_t dim_solver(2*m_ite_max_attempt+1);

    H5std_string GROUP_TEMP;
    hsize_t Max_storage_temp;

    if (solved){
        GROUP_TEMP = "solved";
        Max_storage_temp = m_max_storage_sol;
    }
    else {
        GROUP_TEMP = "unsolved";
        Max_storage_temp = m_max_storage_unsol;
    }

    /*
     * Try block to detect exceptions raised by any of the calls inside it
     */
    try{
        /*
         * Turn off the auto-printing when failure occurs so that we can
         * handle the errors appropriately
         */
        Exception::dontPrint();
        /*
         * Create or Open a file.
         */
        H5File* file;
        try {
            file = new H5File( FILE_NAME, H5F_ACC_RDWR );
        } catch (...) {
            file = new H5File( FILE_NAME, H5F_ACC_TRUNC );

            Group* M_solved = new Group(file->createGroup(GROUP_NAME_I));
            Group(M_solved->createGroup(GROUP_NAME1));
            Group(M_solved->createGroup(GROUP_NAME2));
            Group(M_solved->createGroup(GROUP_NAME3));

            Group* M_unsolved = new Group(file->createGroup(GROUP_NAME_II));
            Group(M_unsolved->createGroup(GROUP_NAME1));
            Group(M_unsolved->createGroup(GROUP_NAME2));
            Group(M_unsolved->createGroup(GROUP_NAME3));

            hsize_t dim_LM[1] = {1};
            DataSpace space_LM( 1, dim_LM );
            DataSet(M_solved->createDataSet( Last_Memb, PredType::NATIVE_INT, space_LM ));
            DataSet(M_unsolved->createDataSet( Last_Memb, PredType::NATIVE_INT, space_LM ));

            hsize_t dim_idx_solver[2] = {1, dim_solver};
            hsize_t maxdims[2] = {H5S_UNLIMITED, dim_solver}; // unlimited dataspace
            DataSpace space_solver( 2, dim_idx_solver, maxdims );
            DSetCreatPropList prop; // Modify dataset creation property to enable chunking
            hsize_t chunk_dims[2] = {1, dim_solver}; // with extendible dataset we cannot use contiguous but chunked dataset
            prop.setChunk(2, chunk_dims);
            DataSet(M_solved->createDataSet( LCP_error, PredType::NATIVE_DOUBLE, space_solver, prop ));
            DataSet(M_unsolved->createDataSet( LCP_error, PredType::NATIVE_DOUBLE, space_solver, prop ));

            delete M_unsolved;
            delete M_solved;
        }
        
        /*
         * Recovering the number of LCP
         */
        Group *Matrix_G, *Vector_G, *Z_G, *Root;

        Group* M_solved = new Group(file->openGroup(GROUP_NAME_I));
        Group* MS = new Group(M_solved->openGroup(GROUP_NAME1));
        hsize_t nb_lcp_sol = MS->getNumObjs();

        Group* M_unsolved = new Group(file->openGroup(GROUP_NAME_II));
        Group* MU = new Group(M_unsolved->openGroup(GROUP_NAME1));
        hsize_t nb_lcp_unsol = MU->getNumObjs();

        hsize_t nb_lcp = nb_lcp_sol + nb_lcp_unsol;

        delete M_unsolved;
        delete M_solved;
        delete MS;
        delete MU;

        /*
         * Checking if the total capacity of storage is reached
         */
        if (int(nb_lcp_sol) > m_max_storage_sol && int(nb_lcp_unsol) > m_max_storage_unsol){
            std::cout << "the maximum storage (" << m_max_storage_sol+m_max_storage_unsol << ") is reached.\n";
            delete file;
            return true;
        }

        Root = new Group(file->openGroup(GROUP_TEMP));
        Matrix_G = new Group(Root->openGroup(GROUP_NAME1));
        Vector_G = new Group(Root->openGroup(GROUP_NAME2));
        Z_G = new Group(Root->openGroup(GROUP_NAME3));

        hsize_t nb_lcp_temp = Matrix_G->getNumObjs();

        if (nb_lcp_temp > Max_storage_temp){
            delete Matrix_G;
            delete Vector_G;
            delete Z_G;
            delete Root;
            delete file;
            return false;
        }

        bool G_exist = 1;
        if (nb_lcp_temp==0) {
            G_exist=0;
            // Initialisation: we fulfill the first line before to extend the dataset
            DataSet* dataset_solver = new DataSet(Root->openDataSet( LCP_error ));
        
            double idx_solv[dim_solver];
            for (hsize_t i=0; i<dim_solver; ++i){
                if (int(i)<2*m_ite_max_attempt) { idx_solv[i] = stats_vec_lcp[i]; }
                else { idx_solv[i] = double(w_fail); }
            }

            dataset_solver->write(idx_solv, PredType::NATIVE_DOUBLE);
        }

        /*
         * Comparison to the previous LCP failure (to prevent similar LCP)
         */
        bool isnt_same_LCP = 1;

        if (G_exist) {
            int last_lcp[1];
            DataSet* dataset_LM = new DataSet(Root->openDataSet( Last_Memb ));
            dataset_LM->read( last_lcp, PredType::NATIVE_INT );
            delete dataset_LM;

            const H5std_string name_data_pre = std::to_string(last_lcp[0]);

            DataSet* dataset_pre = new DataSet(Matrix_G->openDataSet( name_data_pre ));

            DataSpace fspace1 = dataset_pre->getSpace();
            std::size_t dim_out = std::sqrt( fspace1.getSelectNpoints() );
   
            if (dim_out==lcp.dim){
                double data_out[dim_out][dim_out];

                dataset_pre->read( data_out, PredType::NATIVE_DOUBLE );

                /*
                 * Check if matrix already exists? (A large number of attempt to solve LCP)
                 */
                const int dim_D = 3*lcp.dim/4; // size of the Delassus matrix
                Eigen::MatrixXd Diff( dim_D , dim_D );
                for (int i=0; i<dim_D; ++i){
                    for (int j=0; j<dim_D; ++j){
                        const double val_rel = std::min( std::abs(lcp.M(i,j)) , std::abs(data_out[i][j]) );
                        double div = val_rel;
                        if (val_rel==0) {div = 1.0;}
                        const double val_rel_a = (lcp.M(i,j) - data_out[i][j])/div;

                        Diff(i,j) = std::max( std::abs( val_rel_a ) , 0.0);
                    }
                }
                isnt_same_LCP = Diff.norm() > 1e-7;
            } 

            delete dataset_pre;
        }

        /*
         * Create dataspace for the dataset in the file.
         */
        const H5std_string name_matrix = std::to_string(nb_lcp+1);
    
        if (isnt_same_LCP) {
            int dim_M = lcp.dim;
            hsize_t dim_space_M[2];
            dim_space_M[0] = dim_space_M[1] = dim_M;
            DataSpace fspace_M( 2, dim_space_M );
            
            /*
             * Create dataset and write it into the file.
             */
            DataSet* dataset_M = new DataSet(Matrix_G->createDataSet(name_matrix
                , PredType::NATIVE_DOUBLE, fspace_M));

            /*
             * Conversion Eigen -> DOUBLE
             */
            double lcp_M[dim_M][dim_M];
            for (int i=0; i<dim_M; ++i){
                for (int j=0; j<dim_M; ++j){
                    lcp_M[i][j] = lcp.M(i,j);
                }
            }

            dataset_M->write(lcp_M, PredType::NATIVE_DOUBLE);

            /*
             * Close the dataset
             */        
            delete dataset_M;

            /*
             * Save the name of the last written matrix 
             */
            DataSet* dataset_LM = new DataSet(Root->openDataSet( Last_Memb ));
            const int nb_LM[1] = {static_cast<int>(nb_lcp+1)};
            dataset_LM->write(nb_LM, PredType::NATIVE_INT);
            delete dataset_LM;

            /*
             * fulfill extendible dataset
             */
            if (G_exist){
                /*
                 * Save information on solvers with extendible dataset
                 */
                DataSet* dataset_solver = new DataSet(Root->openDataSet( LCP_error ));
                DataSpace space_solver = dataset_solver->getSpace();
                hsize_t dim_curr_s[2]; // dimension of the dataset
                space_solver.getSimpleExtentDims( dim_curr_s, NULL); // retrieves the current dimensions 
                hsize_t ext_size_s[2] = { dim_curr_s[0]+1, dim_curr_s[1]}; 
                dataset_solver->extend( ext_size_s ); // extension with one new line 
      
                DataSpace fspace_s = dataset_solver->getSpace();
                hsize_t dim_s[2] = {1,dim_solver}; 
                hsize_t offset_s[2] = {dim_curr_s[0], 0};
                fspace_s.selectHyperslab( H5S_SELECT_SET, dim_s, offset_s); // selection of the hyperslab
                DataSpace mspace_s( 2, dim_s );
                
                double idx_solv[dim_solver];
                for (hsize_t i=0; i<dim_solver; ++i){
                    if (int(i)<2*m_ite_max_attempt) { idx_solv[i] = stats_vec_lcp[i]; }
                    else { idx_solv[i] = double(w_fail); }
                }
                
                dataset_solver->write(idx_solv, PredType::NATIVE_DOUBLE, mspace_s, fspace_s); // write in the hyperslab
                
                delete dataset_solver;
            }
            /*-----------------------------------------------------------------------------------------
             * new dataset for relative velocities
             *---------------------------------------------------------------------------------------*/
            const H5std_string name_vector = name_matrix;
            /*
             * Create dataspace for the dataset in the file.
             */

            hsize_t dim_space_q[1];
            dim_space_q[0] = dim_M;
            DataSpace fspace_q( 1, dim_space_q );
            /*
             * Create dataset and write it into the file.
             */
            DataSet* dataset_q = new DataSet(Vector_G->createDataSet(name_vector
                , PredType::NATIVE_DOUBLE, fspace_q));

            /*
             * Conversion Eigen -> DOUBLE
             */
            double lcp_q[dim_M];
            for (int i=0; i<dim_M; ++i){
                lcp_q[i] = lcp.q(i);
            }

            dataset_q->write(lcp_q, PredType::NATIVE_DOUBLE);

            /*
             * Close the dataset and the file.
             */        
            delete dataset_q;

            /*-----------------------------------------------------------------------------------------
             * new dataset for z solution of LCP(q,M)
             *---------------------------------------------------------------------------------------*/
            const H5std_string name_z = name_matrix;
            /*
             * Create dataspace for the dataset in the file.
             */

            hsize_t dim_z[1];
            dim_z[0] = dim_M;
            DataSpace fspace_z( 1, dim_z );
            /*
             * Create dataset and write it into the file.
             */
            DataSet* dataset_z = new DataSet(Z_G->createDataSet(name_z
                , PredType::NATIVE_DOUBLE, fspace_z));

            /*
             * Conversion Eigen -> DOUBLE
             */
            double lcp_z[dim_M];
            for (int i=0; i<dim_M; ++i){
                lcp_z[i] = lcp.z(i);
            }
            dataset_z->write(lcp_z, PredType::NATIVE_DOUBLE);

            /*
             * Close the dataset and the file.
             */        
            delete dataset_z;

        }
        delete Matrix_G;
        delete Vector_G;
        delete Z_G;
        delete Root;
        delete file;
   }  // end of try block
   // catch failure caused by the H5File operations
   catch( FileIException error )
   {
    error.printErrorStack();
   }
   // catch failure caused by the DataSet operations
   catch( DataSetIException error )
   {
    error.printErrorStack();
   }
   // catch failure caused by the DataSpace operations
   catch( DataSpaceIException error )
   {
    error.printErrorStack();
   }
   return false;
}

}}} // namespace floe::lcp::solver


#endif // OPE_LCP_SOLVER_HPP