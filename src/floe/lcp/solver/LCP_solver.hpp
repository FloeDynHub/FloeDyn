/*!
 * \file lcp/solver/LCP_solver.hpp
 * \brief LCP solver
 * \author Quentin Jouet
 *//*!
 * \file lcp/solver/LCP_solver.hpp
 * \brief LCP solver
 * \author Quentin Jouet
 */

#ifndef OPE_LCP_SOLVER_HPP
#define OPE_LCP_SOLVER_HPP
// #include <mpi.h>
#include "floe/lcp/solver/LCP_solver.h"

#include "floe/lcp/solver/lexicolemke.hpp"
 // #include "floe/lcp/solver/lexicolemke_eigen.hpp"
#include "floe/lcp/solver/lemke_eigen.hpp"
// #include "floe/lcp/solver/newton_min.hpp"
// #include "floe/lcp/solver/siconos.hpp"
#include "floe/lcp/builder/graph_to_lcp.hpp"
#include "floe/collision/contact_graph.hpp" // boost::edges (todo : should be hidden)
#include <algorithm>
// #include <siconos/lcp_cst.h>
#include <chrono>
// #include <Eigen/SVD>
// #include <JacobiSVD.h>

// saving matrix when lcp solver failed for further analysing
#include "H5Cpp.h"
using namespace H5;

namespace floe { namespace lcp { namespace solver
{


template<typename T>
bool LCPSolver<T>::solve( lcp_type& lcp ) {
        using namespace floe::lcp::solver;
        return (lemke(lcp) || lexicolemke(lcp));
    }


template<typename T>
template<typename TContactGraph>
std::array<vector<typename LCPSolver<T>::real_type>, 2>
LCPSolver<T>::solve( TContactGraph& graph, bool& success, int lcp_failed_stats[] ) {
    using namespace floe::lcp::solver;

    floe::lcp::builder::GraphLCP<real_type, decltype(graph)> graph_lcp( graph );
    auto lcp_orig = graph_lcp.getLCP();
    auto lcp = lcp_orig;
    std::vector<int> lcp_test_list{1,2,3,3,3,4,4,4}; //lcp_test_list{4,4,4}; // best config: lcp_test_list{1,2,3,3,3,4,4,4};

    int solver_used, test_used{0};
    static bool is_full_storage = false;
    real_type Err;
    real_type ECc;
    bool vitrelnormtest; 
    int w_fail;

    // %%%%%%%%%%%%%%%%%%%%%%%%
    // % phase de compression %
    // %%%%%%%%%%%%%%%%%%%%%%%%

    bool solved{0};
    vector<real_type> Solc(3 * graph_lcp.nb_floes), floe_impulses(graph_lcp.nb_floes, 0);

    decltype(lcp.z) best_z;
    real_type best_Err = std::numeric_limits<real_type>::max();
    bool solver_success{0};

    // for (int i = 0; i < 60000; ++i) lemke(lcp);

    // int nb_solvers{5};
    for (auto test_idx : lcp_test_list)
    {
        ++test_used;

        for (int i = 0; i < m_nb_solvers; ++i)
        {
            // int mpi_rk;
            // MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rk );
            
            solver_used = i+1;
 
            if (test_used<=5 && i==1) {break;}// for best config: do not take into account lemke before case 4.
            solver_success = run_solver(lcp, i);

            // if (!solver_success)
            //     continue;

            if (std::any_of(lcp.z.begin(), lcp.z.end(), is_nan<double>))
                { std::cout << "***NAN***"; continue; }

            lcp_orig.z = lcp.z;

             // As-t-on une meilleure solution ?
            Err = LCP_error(lcp_orig);

            if (Err < best_Err || best_Err == std::numeric_limits<real_type>::max())
            {
                best_z = lcp.z;
                best_Err = Err;
            }

             // the best solution is always tested
            lcp_orig.z = best_z;
            Err = best_Err;

             // Solution correspondante
            Solc = calcSolc(graph_lcp, lcp_orig);

            if (std::any_of(Solc.begin(), Solc.end(), is_nan<double>))
                continue;

             // Energie cinetique, Erreur LCP & Vit rel Normale :
            ECc = calcEc(Solc, graph_lcp.M, graph_lcp.W);

            vitrelnormtest = VRelNtest(prod(trans(graph_lcp.J), Solc), graph);

            solved = LCPtest(test_idx,ECc,1,Err,vitrelnormtest);

            // if (lcp.dim > 200){ std::cout << "lcp-rk-" << mpi_rk << "-dim-" << lcp.dim << "-tidx-" << test_idx << "-success-" << solved << ", " << std::flush; }
            if (solved) {
                // m_solver_stats(test_idx-1, i) += 1; // test
                // if (!solver_success) std::cout << "SOLVER_LIE(" << i << ")"; // test
                break; }
        }
        if (solved) {
            if (!is_full_storage){
                w_fail = which_failure(Err,ECc,vitrelnormtest,lcp_test_list[test_used-1]);
                is_full_storage = saving_LCP_in_hdf5(lcp_orig, solved, test_used, solver_used, Err, w_fail );
            }

            break;
        }
        // lcp = random_perturbation(lcp, 1e-10);
        // auto lcp = lcp_orig;
        // auto Err = LCP_error(lcp_orig);
        // auto vitrelnormtest = VRelNtest(prod(trans(graph_lcp.J), Solc), graph);
        // auto ECc = calcEc(Solc, graph_lcp.M, graph_lcp.W);

        // if (test_idx==4){
        //     random_perturbation2(lcp, 5*1e-11); //// influence of the perturbations
        // }
        random_perturbation2(lcp, 5*1e-11); //// influence of the perturbations

    }
    
    if (!solved) {
        //Matt
        if (!is_full_storage){
            w_fail = which_failure(Err,ECc,vitrelnormtest,lcp_test_list[test_used-1]);
            is_full_storage = saving_LCP_in_hdf5(lcp_orig, solved, test_used, solver_used, Err, w_fail );
        }

        lcp_failed_stats[0] += 1;
        // EndMatt

        success = 0;
        return {{graph_lcp.W, floe_impulses}};
    }

    // %%%%%%%%%%%%%%%%%%%%%%%%%%
    // % phase de decompression %
    // %%%%%%%%%%%%%%%%%%%%%%%%%%
    vector<real_type> Sold(3 * graph_lcp.nb_floes);
    if ( epsilon != 0 )
    {
        lcp_type lcp_d_orig = lcp = graph_lcp.getLCP_d(lcp_orig, Solc, epsilon);
        real_type born_sup = graph_lcp.born_sup_d(lcp_orig, epsilon);
        // if (born_sup > 1) std::cout << "*****BORN SUP***** " << born_sup << " ";

        solved = 0;
        best_Err = std::numeric_limits<real_type>::max();

        for (auto test_idx : lcp_test_list)
        {
            for (int i = 0; i < m_nb_solvers; ++i)
            {
                solver_success = run_solver(lcp, i);

                // if (!solver_success)
                //     continue;

                if (std::any_of(lcp.z.begin(), lcp.z.end(), is_nan<double>))
                    continue;

                lcp_d_orig.z = lcp.z;

                 // As-t-on une meilleure solution ?
                // auto Err = LCP_error(lcp_d_orig);
                Err = LCP_error(lcp_d_orig);
                if (Err < best_Err || best_Err == std::numeric_limits<real_type>::max())
                {
                    best_z = lcp.z;
                    best_Err = Err;
                }

                 // On teste toujours la meilleure solution trouvée
                lcp_d_orig.z = best_z;
                Err = best_Err;
                
                 // Solution correspondante
                Sold = calcSold(graph_lcp, lcp_orig, lcp_d_orig, Solc);

                if (std::any_of(Sold.begin(), Sold.end(), is_nan<double>))
                    continue; // std::cout << "******NAN******";

                 // Energie cinetique, Erreur LCP & Vit rel Normale :
                auto ECd = calcEc(Sold, graph_lcp.M, graph_lcp.W);
                // Err = LCP_error(lcp_d_orig);
                // auto vitrelnormtest = VRelNtest(prod(trans(graph_lcp.J), Sold), graph);
                vitrelnormtest = VRelNtest(prod(trans(graph_lcp.J), Sold), graph);
                solved = LCPtest(test_idx, ECd, 1 + born_sup, Err, vitrelnormtest);

                //solution alternative pour la conservation de l'EC
                if (solved)
                {
                    if (ECd > 1)
                    {
                        lcp_failed_stats[2] += 1;
                        Sold = (1 + epsilon) * Solc - epsilon * graph_lcp.W; // return this instead of Sold
                        floe_impulses = graph_lcp.impulse_vector(lcp_orig, epsilon);
                    } else {
                        // Impulse calculation
                        floe_impulses = graph_lcp.impulse_vector(lcp_orig, lcp_d_orig, epsilon);
                    }
                    break;
                }
            }
            if (solved) break;
            // lcp = random_perturbation(lcp, 1e-10);
            random_perturbation2(lcp, 5*1e-11);
        }
        if (!solved)
        {
            //Matt
            lcp_failed_stats[1] += 1;
            //EndMatt
            success = 0;
            floe_impulses = graph_lcp.impulse_vector(lcp_orig, epsilon);
            return {{(1 + epsilon) * Solc - epsilon * graph_lcp.W, floe_impulses}};
        }
    } else {
        success = 1;
        floe_impulses = graph_lcp.impulse_vector(lcp_orig);
        return {{Solc, floe_impulses}};
    }

    success = 1;
    return {{Sold, floe_impulses}};

}

template<typename T>
bool LCPSolver<T>::run_solver(lcp_type& lcp, int id){
    nb_solver_run++;
    auto t_start = std::chrono::high_resolution_clock::now(); // test
    switch (id)
    {
        // case 0:
        //     return lemke(lcp);
        // case 1:
        //     return lexicolemke(lcp);
        // case 2:
        //     // lcp.z = zero_vector<T>(lcp.dim);
        //     return true;

        // case 0:
        //     // lcp.z = zero_vector<T>(lcp.dim);
        //     return true;
        // case 1:
        //     return lemke(lcp);
        // case 2:
        //     return lexicolemke(lcp);

        case 0:
            lexicolemke(lcp);
            break;
            // return true;
        case 1:
            lemke(lcp);
            break;
            // return true;

        // case 2:
        //     // lcp.z = zero_vector<T>(lcp.dim);
        //     return true;
        // case 2:
        //     {
        //     // auto lcp_copy = lcp; // newton_min seems to modify lcp (TODO test again)
        //     // bool success = newton_min(lcp_copy);
        //     // lcp.z = lcp_copy.z;
        //     // return success;
        //         // return true;
        //     return newton_min(lcp);
        //     }
        // case 3:
        //     return siconos(lcp, SICONOS_LCP_RPGS);
        // case 4:
        //     return siconos(lcp, SICONOS_LCP_AVI_CAOFERRIS);
        // case 5:
        //     return siconos(lcp, SICONOS_LCP_NSQP);
        // case 6:
        //     return siconos(lcp, SICONOS_LCP_PIVOT);
        // case 6:
        //     return siconos(lcp, SICONOS_LCP_CPG);
        // case 7:
        //     return siconos(lcp, SICONOS_LCP_CPG);
        // case 8:
        //     return siconos(lcp, SICONOS_LCP_CPG);
    }
    auto t_end = std::chrono::high_resolution_clock::now(); // test
    auto call_time = std::chrono::duration<double, std::milli>(t_end-t_start).count(); // test
    chrono_solver += call_time;
    max_chrono_solver = std::max(max_chrono_solver, call_time);
    return 1;
}


template<typename T>
bool LCPSolver<T>::LCPtest(int compt, real_type EC, real_type born_EC, real_type Err, bool VRelNtest ){
    if (compt == 1)
    {
        if (EC > born_EC*(1+1e-4) || Err > 5*1e-11 || !VRelNtest)
            return false;
    }
    else if (compt == 2)
    {
        if (EC > born_EC * (1+1e-4) || Err > 1e-8 || !VRelNtest)
            return false;
    }
    else if (compt == 3)
    {
        if (EC > born_EC*(1+1e-2) || Err > 1e-2 || !VRelNtest)
            return false;
    }
    else if (compt == 4)
    {
        if (EC > born_EC*(1+1e-2) || !VRelNtest) // no error test
            return false;
    }
    return true;
}


template<typename T>
template<typename Tmat, typename Tvect>
typename LCPSolver<T>::real_type 
LCPSolver<T>::calcEc(const Tvect& S, const Tmat& M, const Tvect& w){
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
typename LCPSolver<T>::lcp_type
LCPSolver<T>::random_perturbation(lcp_type& lcp, real_type max){
    // version randomly modifying all non-zeros values
    for (auto it1 = lcp.A.begin1(); it1 != lcp.A.end1(); ++it1) {
        for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
            if (*it2 != 0)
                *it2 += random_real(max);
        }
    }
    /* // Matlab version (divmat), modifying only one value in lcp.A and one in lcp.q (worse results)
    auto m = lcp.dim / 4;
    std::size_t n1 = std::rand() % (3 * m), n2 = std::rand() % (3 * m); // Iheart variable (=3*m) in Matlab code
    lcp.A(n1,n2) *= (1 + random_real(max));
    lcp.q(n1) *= (1 + random_real(max));
    */
    return lcp; 
}

template<typename T>
void LCPSolver<T>::random_perturbation2(lcp_type& lcp, real_type max){
    // version randomly modifying all non-zeros values
    for (auto it1 = lcp.A.begin1(); it1 != lcp.A.end1(); ++it1) {
        for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
            if (*it2 != 0)
                *it2 += random_real(max);
                // if (lcp.dim<9){std::cout << "random term:" << random_real(max);}
        }
    }
    /* // Matlab version (divmat), modifying only one value in lcp.A and one in lcp.q (worse results)
    auto m = lcp.dim / 4;
    std::size_t n1 = std::rand() % (3 * m), n2 = std::rand() % (3 * m); // Iheart variable (=3*m) in Matlab code
    lcp.A(n1,n2) *= (1 + random_real(max));
    lcp.q(n1) *= (1 + random_real(max));
    */
}

template<typename T>
template<typename TContactGraph>
bool LCPSolver<T>::VRelNtest(const vector<real_type>& V, const TContactGraph& graph){
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
T LCPSolver<T>::random_real(T max)
{
    return m_uniform_distribution(m_random_generator) * max;
}


template<typename T>
bool LCPSolver<T>::saving_LCP_in_hdf5(lcp_type& lcp, bool solved, int test_idx, int solver_used, real_type lcp_err, int w_fail){

    const H5std_string FILE_NAME("/Users/matthiasrabatel/Travail/outputs_mycode/matrix.h5");
    const H5std_string GROUP_NAME_I( "solved" ); // root group
    const H5std_string GROUP_NAME_II( "unsolved" ); // root group
    const H5std_string GROUP_NAME1( "M" ); 
    const H5std_string GROUP_NAME2( "q" );
    const H5std_string GROUP_NAME3( "z" );
    const H5std_string LCP_error( "LCP error" );
    const H5std_string Last_Memb( "Last LCP" ); // to prevent similar LCP
    const H5std_string Contact_Graph_Info( "Contact Graph Info" ); // to store information on the contact graph and the "while loop"
    const H5std_string Idx_solver( "Which solver" ); // Information on which solver, 
    // how many random perturbations are used before to compute solution, the index in the h5 file and
    // the source of the LCP error (see which_failure)
    const hsize_t Max_storage_sol = 15000;
    const hsize_t Max_storage_unsol = 15000;


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

            hsize_t dim_idx_solver[2] = {1, 4};
            hsize_t maxdims[2] = {H5S_UNLIMITED, 4}; // unlimited dataspace
            DataSpace space_solver( 2, dim_idx_solver, maxdims );
            DSetCreatPropList prop; // Modify dataset creation property to enable chunking
            hsize_t chunk_dims[2] = {1, 4}; // with extendible dataset we cannot use contiguous but chunked dataset
            prop.setChunk(2, chunk_dims);
            DataSet(M_solved->createDataSet( Idx_solver, PredType::NATIVE_INT, space_solver, prop ));
            DataSet(M_unsolved->createDataSet( Idx_solver, PredType::NATIVE_INT, space_solver, prop ));

            hsize_t maxdims_le[1] = {H5S_UNLIMITED};
            DataSpace space_LE( 1, dim_LM, maxdims_le );
            DSetCreatPropList prop_le; // Modify dataset creation property to enable chunking
            hsize_t chunk_dims_le[1] = {1}; // with extendible dataset we cannot use contiguous but chunked dataset
            prop_le.setChunk(1, chunk_dims_le);
            DataSet(M_solved->createDataSet( LCP_error, PredType::NATIVE_DOUBLE, space_LE, prop_le ));
            DataSet(M_unsolved->createDataSet( LCP_error, PredType::NATIVE_DOUBLE, space_LE, prop_le ));

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
            DataSet* dataset_solver = new DataSet(Root->openDataSet( Idx_solver ));
            int idx_solv[4] = { test_idx , solver_used,  static_cast<int>(nb_lcp+1), w_fail };
            dataset_solver->write(idx_solv, PredType::NATIVE_INT);

            DataSet* dataset_LE = new DataSet(Root->openDataSet( LCP_error ));
            dataset_LE->write(&lcp_err, PredType::NATIVE_DOUBLE);
        }

        /*
         * Comparison to the previous LCP failure (to prevent similar LCP)
         */
        bool isnt_same_LCP = 1;

        if (G_exist) {
            int last_lcp[1];
            DataSet* dataset_LM = new DataSet(Root->openDataSet( Last_Memb ));
            dataset_LM->read( last_lcp, PredType::NATIVE_INT );

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
                        const double val_rel = std::min( std::abs(lcp.A(i,j)) , std::abs(data_out[i][j]) );
                        double div = val_rel;
                        if (val_rel==0) {div = 1.0;}
                        const double val_rel_a = (lcp.A(i,j) - data_out[i][j])/div;

                        Diff(i,j) = std::max( std::abs( val_rel_a ) , 0.0);
                    }
                }
                isnt_same_LCP = Diff.norm() > 1e-7;
            } 
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
                    lcp_M[i][j] = lcp.A(i,j);
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
                DataSet* dataset_solver = new DataSet(Root->openDataSet( Idx_solver ));
                DataSpace space_solver = dataset_solver->getSpace();
                hsize_t dim_curr_s[2]; // dimension of the dataset
                space_solver.getSimpleExtentDims( dim_curr_s, NULL); // retrieves the current dimensions 
                hsize_t ext_size_s[2] = { dim_curr_s[0]+1, dim_curr_s[1]}; 
                dataset_solver->extend( ext_size_s ); // extension with one new line 
      
                DataSpace fspace_s = dataset_solver->getSpace();
                hsize_t dim_s[2] = {1,4}; 
                hsize_t offset_s[2] = {dim_curr_s[0], 0};
                fspace_s.selectHyperslab( H5S_SELECT_SET, dim_s, offset_s); // selection of the hyperslab
                DataSpace mspace_s( 2, dim_s );
                int idx_solv[4] = { test_idx , solver_used, static_cast<int>(nb_lcp+1), w_fail };
                dataset_solver->write(idx_solv, PredType::NATIVE_INT, mspace_s, fspace_s); // write in the hyperslab
                
                delete dataset_solver;

                /*
                 * Save information on LCP error with extendible dataset
                 */                
                DataSet* dataset_LE = new DataSet(Root->openDataSet( LCP_error ));
                DataSpace space_LE = dataset_LE->getSpace();
                hsize_t dim_curr_le[1]; // dimension of the dataset
                space_LE.getSimpleExtentDims( dim_curr_le, NULL); // retrieves the current dimensions 
                hsize_t ext_size_le[1] = { dim_curr_le[0]+1}; 
                dataset_LE->extend( ext_size_le ); // extension with one new line 
      
                DataSpace fspace_le = dataset_LE->getSpace();
                hsize_t dim_le[1] = {1}; 
                hsize_t offset_le[1] = {dim_curr_le[0]};
                fspace_le.selectHyperslab( H5S_SELECT_SET, dim_le, offset_le); // selection of the hyperslab
                DataSpace mspace_le( 1, dim_le );
                dataset_LE->write(&lcp_err, PredType::NATIVE_DOUBLE, mspace_le, fspace_le); // write in the hyperslab

                delete dataset_LE;
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


template<typename T>
int LCPSolver<T>::which_failure(real_type Err, real_type EC, bool rel_n_vel, int compt){

    int source=0;
    real_type err_limit{0}, ec_limit{0};
    if (compt==1) {err_limit=5*1e-11; ec_limit=1+1e-4;}
    else if (compt==2) {err_limit=1e-8; ec_limit=1+1e-4;}
    else if (compt==3) {err_limit=1e-2; ec_limit=1+1e-2;}
    else if (compt==4) {err_limit=std::numeric_limits<real_type>::max(); ec_limit=1+1e-2;}

    if (Err>err_limit){source += 100;}
    if (EC>ec_limit) {source += 20;}
    if (!rel_n_vel){source += 3;}

    return source;
}

}}} // namespace floe::lcp::solver


#endif // OPE_LCP_SOLVER_HPP