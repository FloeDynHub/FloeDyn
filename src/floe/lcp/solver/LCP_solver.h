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

using namespace H5;

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

    LCPSolver(real_type epsilon) : epsilon{epsilon}, tolerance{1e-7}, coef_perturb{1e-9}, 
        ite_max_attempt{5} {}

    template<typename TContactGraph>
    std::array<vector<real_type>, 2> solve( TContactGraph& graph, bool& success, int lcp_failed_stats[] );

    double chrono_solver{0.0}; // test perf
    double max_chrono_solver{0.0}; // test perf

protected:
    typedef boost::numeric::ublas::matrix<real_type> array_type;
    typedef boost::numeric::ublas::matrix<real_type> vector_type;

    real_type epsilon;          //!< energy restitution coeff
    real_type tolerance;        // tolerance consented for the LCP error
    real_type coef_perturb;     // multiplier coefficient for the size of the perturbations (1e-9 seems to be an optimal value). 
                                // One could be taking from 1e-7 to 1e-10.
    int       ite_max_attempt;  // integer for the number of perturbations (5 seems to be a good compromise between 
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

////////////////////////////////////////////////////////////
//! Random small perturbation of LCP: see #5 and matlab file to change the random perturbation routines. 
// Keeping only the case 2/ and 3/ (=reduction_via_perturbation) to only use on IterLemke and Lexico 
// with a coef setted to 1e-8.
////////////////////////////////////////////////////////////
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

/*  \fn saving_LCP_in_hdf5(lcp_type lcp, bool solved, int count_attempt, int count_RP, 
 *      int count_SR, int count_SR_failed, const int perturb_used, bool use_lexico_ordering, 
 *      real_type lcp_err, int w_fail)
 *  
 *  \brief  Saving M and q from dealt with LCP(M,q) (solved and unsolved) for further analysis
 *          Return a boolean to prevent the maximum capacity to store (ex: 50 000 LCP ~ 250 Mo).
 *
 *  \remark The storage used hdf5 file formulation and consists in the following kind of table:
 *          |     1     |     2      |     3      |   4   |          5         |   6    |      8      |
 *          |:---------:|:----------:|:----------:|:-----:|:------------------:|:------:|:-----------:|
 *          | LCP error | nb attempt | nb perturb | nb SR | nb adj cone failed | lexico | idx failure |
 *
 *          with: \e nb for number, \e SR for secondary ray, \e adj \e cone \e failed for the failure of the 
 *          method consisting in going through an adjacent cone, \e lexico is true if the lexicographic 
 *          ordering is used during at least one Lemke's algorithm, \e idx \e perturb for the index of the
 *          matrix perturbation (see LCPSolver::matrix_perturbation(const std::size_t dim , matrix &M, const double alpha, const int Idx_perturb))
 *          and \e idx \e failure for the source of the LCP error (see LCPSolver::which_failure( vector<real_type> Err, bool Is_pos_rel_norm_vel )).
 */
template<typename T>
bool saving_LCP_in_hdf5(floe::lcp::LCP<T> lcp, bool solved, int count_attempt, int count_RP, 
    int count_SR, int count_SR_failed, bool use_lexico_ordering, T lcp_err, int w_fail) 
{
    const H5std_string FILE_NAME("/Users/matthiasrabatel/Travail/outputs_mycode/matrix.h5");
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
    const hsize_t Max_storage_sol = 15000;
    const hsize_t Max_storage_unsol = 15000;

    const hsize_t dim_solver(7);

    H5std_string GROUP_TEMP;
    hsize_t Max_storage_temp;

    if (solved){
        GROUP_TEMP = "solved";
        Max_storage_temp = Max_storage_sol;
    }
    else {
        GROUP_TEMP = "unsolved";
        Max_storage_temp = Max_storage_unsol;
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

            // hsize_t maxdims_le[1] = {H5S_UNLIMITED};
            // DataSpace space_LE( 1, dim_LM, maxdims_le );
            // DSetCreatPropList prop_le; // Modify dataset creation property to enable chunking
            // hsize_t chunk_dims_le[1] = {1}; // with extendible dataset we cannot use contiguous but chunked dataset
            // prop_le.setChunk(1, chunk_dims_le);
            // DataSet(M_solved->createDataSet( LCP_error, PredType::NATIVE_DOUBLE, space_LE, prop_le ));
            // DataSet(M_unsolved->createDataSet( LCP_error, PredType::NATIVE_DOUBLE, space_LE, prop_le ));

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
        if (nb_lcp_sol > Max_storage_sol && nb_lcp_unsol > Max_storage_unsol){
            std::cout << "the maximum storage (" << Max_storage_sol+Max_storage_unsol << ") is reached.\n";
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
            double idx_solv[dim_solver] = { lcp_err, double(count_attempt), double(count_RP), 
                double(count_SR), double(count_SR_failed), double(use_lexico_ordering), 
                double(w_fail) };
            dataset_solver->write(idx_solv, PredType::NATIVE_DOUBLE);

            // DataSet* dataset_LE = new DataSet(Root->openDataSet( LCP_error ));
            
            // dataset_LE->write( &lcp_err, PredType::NATIVE_DOUBLE );
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
                double idx_solv[dim_solver] = { lcp_err, double(count_attempt), double(count_RP), 
                    double(count_SR), double(count_SR_failed), double(use_lexico_ordering), 
                    double(w_fail) };
                // int idx_solv[dim_solver] = { test_idx , solver_used, perturb_used, static_cast<int>(nb_lcp+1), w_fail, 
                //     min_how_is_solved };
                dataset_solver->write(idx_solv, PredType::NATIVE_DOUBLE, mspace_s, fspace_s); // write in the hyperslab
                
                delete dataset_solver;

                // /*
                //  * Save information on LCP error with extendible dataset
                //  */                
                // DataSet* dataset_LE = new DataSet(Root->openDataSet( LCP_error ));
                // DataSpace space_LE = dataset_LE->getSpace();
                // hsize_t dim_curr_le[1]; // dimension of the dataset
                // space_LE.getSimpleExtentDims( dim_curr_le, NULL); // retrieves the current dimensions 
                // hsize_t ext_size_le[1] = { dim_curr_le[0]+1}; 
                // dataset_LE->extend( ext_size_le ); // extension with one new line 
      
                // DataSpace fspace_le = dataset_LE->getSpace();
                // hsize_t dim_le[1] = {1}; 
                // hsize_t offset_le[1] = {dim_curr_le[0]};
                // fspace_le.selectHyperslab( H5S_SELECT_SET, dim_le, offset_le); // selection of the hyperslab
                // DataSpace mspace_le( 1, dim_le );

                // dataset_LE->write( &lcp_err, PredType::NATIVE_DOUBLE, mspace_le, fspace_le); // write in the hyperslab

                // delete dataset_LE;
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
    error.printError();
   }
   // catch failure caused by the DataSet operations
   catch( DataSetIException error )
   {
    error.printError();
   }
   // catch failure caused by the DataSpace operations
   catch( DataSpaceIException error )
   {
    error.printError();
   }
   return false;
}

}}} // namespace floe::lcp::solver


#endif // OPE_LCP_SOLVER_H