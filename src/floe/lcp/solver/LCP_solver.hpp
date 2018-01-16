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
#include <Eigen/SVD>
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
LCPSolver<T>::solve( TContactGraph& graph, bool& success, double lcp_failed_stats[] ) {
    using namespace floe::lcp::solver;

    floe::lcp::builder::GraphLCP<real_type, decltype(graph)> graph_lcp( graph );
    auto lcp_orig = graph_lcp.getLCP();
    auto lcp = lcp_orig;
    std::vector<int> lcp_test_list{1, 2, 3};

    int solver_used;

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
        for (int i = 0; i < m_nb_solvers; ++i)
        {
            // int mpi_rk;
            // MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rk );
            
            solver_used = i;

            solver_success = run_solver(lcp, i);

            // if (!solver_success)
            //     continue;

            if (std::any_of(lcp.z.begin(), lcp.z.end(), is_nan<double>))
                { std::cout << "***NAN***"; continue; }

            lcp_orig.z = lcp.z;

             // As-t-on une meilleure solution ?
            auto Err = LCP_error(lcp_orig);
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
            auto ECc = calcEc(Solc, graph_lcp.M, graph_lcp.W);
            // Err = LCP_error(lcp_orig);
            auto vitrelnormtest = VRelNtest(prod(trans(graph_lcp.J), Solc), graph);
            solved = LCPtest(test_idx,ECc,1,Err,vitrelnormtest);

            // if (lcp_orig.dim < 9){
            //     auto JW = prod(trans(graph_lcp.J), Solc);
            //     real_type minJW = std::numeric_limits<real_type>::max();
            //     for (int j=0; j<graph_lcp.nb_contacts; ++j){
            //         if (JW(j) < minJW){minJW = JW(i);}
            //     }
            //     std::cout << "before random: \n" 
            //         << "A: " << lcp_orig.A << ", \n"
            //         << "q: " << lcp_orig.q << ", \n"
            //         << "z: " << lcp_orig.z << ", \n"
            //         << "min of J^TW: " << minJW <<", \n"
            //         << "VN: " << vitrelnormtest  <<", \n"
            //         << "ECc: " << ECc <<", \n"
            //         << "Err: " << Err <<", \n"
            //         << "test_idx: " << test_idx <<", \n"
            //         << "solver_used: " << solver_used <<", \n"
            //         << "solved: " << solved <<", \n"
            //         << "------------------------------------------ \n";
            // }

            // if (lcp.dim > 200){ std::cout << "lcp-rk-" << mpi_rk << "-dim-" << lcp.dim << "-tidx-" << test_idx << "-success-" << solved << ", " << std::flush; }
            if (solved) {
                m_solver_stats(test_idx-1, i) += 1; // test
                // if (!solver_success) std::cout << "SOLVER_LIE(" << i << ")"; // test
                break; }
        }
        if (solved) {
            // int size_delassus = 3*lcp_orig.dim/4;
            // Eigen::MatrixXd BMB(size_delassus,size_delassus);
            // for (int i=0;i<size_delassus;++i){
            //     for (int j=0;j<size_delassus;++j){
            //         BMB(i,j) = lcp_orig.A(i,j);
            //     }
            // }
            // Eigen::JacobiSVD<Eigen::MatrixXd> svd(BMB);
            // double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
            // if (cond<1e28){
            //     lcp_failed_stats[3] += cond;
            // }
            // else{
            //     lcp_failed_stats[6] += 1;
            // }
            // lcp_failed_stats[3] += min(prod(trans(graph_lcp.J), Solc));

            break;
        }
        // lcp = random_perturbation(lcp, 1e-10);
        // auto lcp = lcp_orig;
        // auto Err = LCP_error(lcp_orig);
        // auto vitrelnormtest = VRelNtest(prod(trans(graph_lcp.J), Solc), graph);
        // auto ECc = calcEc(Solc, graph_lcp.M, graph_lcp.W);

        random_perturbation2(lcp, 5*1e-11);

        // if (lcp_orig.dim < 11){
        //     std::cout << "after random: \n" 
        //             << "A: " << lcp.A << ", \n"
        //             << "q: " << lcp.q << ", \n"
        //             << "z: " << lcp.z << ", \n"
        //             << "J^TW: " << prod(trans(graph_lcp.J), Solc) <<", \n"
        //             << "------------------------------------------ \n";
        // }
    }
    
    if (!solved) {
        //Matt
        saving_matrix_unsolved_LCP(lcp_orig);

        lcp_failed_stats[0] += 1;

        auto Err = LCP_error(lcp_orig);

        auto Err_Coul = LCP_err_Coul(lcp_orig);

        auto ECc = calcEc(Solc, graph_lcp.M, graph_lcp.W);

        auto JW = prod(trans(graph_lcp.J), Solc);
        real_type minJW = std::numeric_limits<real_type>::max();
        for (int j=0; j<graph_lcp.nb_contacts; ++j){
            if (JW(j) < minJW){minJW = JW(j);}
        }

        real_type minECoul = std::numeric_limits<real_type>::max();
        for (int j=0; j<graph_lcp.nb_contacts; ++j){
            if (Err_Coul(j) < minECoul){minECoul = Err_Coul(j);}
        }

        if (ECc>1 && Err <= 1e-2){lcp_failed_stats[3] += 1;
            std::cout << "ECc = " << ECc << ",\n"
                << "Err = " << Err << ", \n";
        }
        if (minJW<0){
            lcp_failed_stats[4] += 1;
            lcp_failed_stats[5] += minJW;
        }
        if (minECoul<0){
            lcp_failed_stats[6] += 1;
            lcp_failed_stats[7] += minECoul;
        }

        lcp_failed_stats[8] += Err;

        // if (lcp_orig.dim < 9 && Err > 10) {
        //     int size_delassus = 3*lcp_orig.dim/4;
        //     Eigen::MatrixXd BMB(size_delassus,size_delassus);
        //     for (int i=0;i<size_delassus;++i){
        //         for (int j=0;j<size_delassus;++j){
        //             BMB(i,j) = lcp_orig.A(i,j);
        //         }
        //     }

        //     std::cout << "A: " << lcp_orig.A << ", \n"
        //         << "q: " << lcp_orig.q << ", \n"
        //         << "z: " << lcp_orig.z << ", \n"
        //         << "min of J^TW: " << minJW <<", \n"
        //         << "ECc: " << ECc <<", \n"
        //         << "Err: " << Err <<", \n"
        //         << "minECoul: " << minECoul <<", \n"
        //         << "------------------------------------------ \n";
        // }
        // int size_delassus = 3*lcp_orig.dim/4;
        // Eigen::MatrixXd BMB(size_delassus,size_delassus);
        // for (int i=0;i<size_delassus;++i){
        //     for (int j=0;j<size_delassus;++j){
        //         BMB(i,j) = lcp_orig.A(i,j);
        //     }
        // }
        // Eigen::JacobiSVD<Eigen::MatrixXd> svd(BMB);
        // double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
        // if (cond<1e28){
        //     lcp_failed_stats[4] += cond;
        // }
        // else{
        //     lcp_failed_stats[7] += 1;
        // }

        // if (Err>1 && lcp.dim < 11){
        //     std::cout << "phase compression: \n" 
        //             << "A: " << lcp_orig.A << ", \n"
        //             << "Cond A: " << cond << ", \n"
        //             << "q: " << lcp_orig.q << ", \n"
        //             << "z: " << lcp_orig.z << ", \n"
        //             << "Err:" << Err << ", \n"
        //             << "------------------------------------------ \n";
        // }
        //EndMatt
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
                auto Err = LCP_error(lcp_d_orig);
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
                auto vitrelnormtest = VRelNtest(prod(trans(graph_lcp.J), Sold), graph);
                solved = LCPtest(test_idx, ECd, 1 + born_sup, Err, vitrelnormtest);

                //solution alternative pour la conservation de l'EC
                if (solved)
                {
                    if (ECd > 1)
                    {
                        lcp_failed_stats[2] += 1;
                        std::cout << "ECD>1 => " << ECd << " ";
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
            // auto Err = LCP_error(lcp_d_orig);
            lcp_failed_stats[1] += 1;

            // int size_delassus = 3*lcp_d_orig.dim/4;
            // Eigen::MatrixXd BMB(size_delassus,size_delassus);
            // for (int i=0;i<size_delassus;++i){
            //     for (int j=0;j<size_delassus;++j){
            //         BMB(i,j) = lcp_d_orig.A(i,j);
            //     }
            // }
            // Eigen::JacobiSVD<Eigen::MatrixXd> svd(BMB);
            // double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
            // if (cond<1e28){
            //     lcp_failed_stats[5] += cond;
            // }
            // else{
            //     lcp_failed_stats[8] += 1;
            // }

            // if (Err>1 && lcp.dim < 11){
            //     std::cout << "phase decompression: \n" 
            //             << "A: " << lcp_d_orig.A << ", \n"
            //             << "Cond A: " << cond << ", \n"
            //             << "q: " << lcp_d_orig.q << ", \n"
            //             << "z: " << lcp_d_orig.z << ", \n"
            //             << "Err:" << Err << ", \n"
            //             << "------------------------------------------ \n";
            // }
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
void LCPSolver<T>::saving_matrix_unsolved_LCP(lcp_type& lcp){
    const H5std_string FILE_NAME( "/Users/matthiasrabatel/Travail/outputs_mycode/matrix.h5" );
    const H5std_string GROUP_NAME1( "Delassus Matrix of unsolved LCP" );
    const H5std_string GROUP_NAME2( "Corresponding relative velocities" );

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
        }
        
        /*
         * Create or Open the groups
         */
        Group* Matrix_G;
        Group* Vector_G;
        try {
            Matrix_G = new Group(file->openGroup(GROUP_NAME1));
            Vector_G = new Group(file->openGroup(GROUP_NAME2));
        } catch (...) {
            /* Create group for floe shapes */
            Matrix_G = new Group(file->createGroup(GROUP_NAME1));
            Vector_G = new Group(file->createGroup(GROUP_NAME2));
        }

        hsize_t Group_size = Matrix_G->getNumObjs();
        const H5std_string name_matrix = std::to_string(Group_size+1);
        /*
         * Create dataspace for the dataset in the file.
         */
        const int dim_M = lcp.dim;
        hsize_t dim_space_M[2];
        dim_space_M[0] = dim_M;
        dim_space_M[1] = dim_M;
        DataSpace fspace_M( 2, dim_space_M );
        /*
         * Create dataset and write it into the file.
         */
        DataSet* dataset_M = new DataSet(Matrix_G->createDataSet(name_matrix
            , PredType::NATIVE_DOUBLE, fspace_M));

        /*
         * Conversion Eigen -> DOUBLE
         */
        double Delassus[dim_M][dim_M];
        for (int i=0; i<dim_M; ++i){
            for (int j=0; j<dim_M; ++j){
                Delassus[i][j] = lcp.A(i,j);
            }
        }

        dataset_M->write(Delassus, PredType::NATIVE_DOUBLE);

        /*
         * Close the dataset and the file.
         */        
        delete dataset_M;

        /*-----------------------------------------------------------------------------------------
         * new dataset for relative velocities
         *---------------------------------------------------------------------------------------*/
        const H5std_string name_vector = std::to_string(Group_size+1);
        /*
         * Create dataspace for the dataset in the file.
         */

        hsize_t dim_space_V[1];
        dim_space_V[0] = dim_M;

        DataSpace fspace_V( 1, dim_space_V );
        /*
         * Create dataset and write it into the file.
         */

        DataSet* dataset_V = new DataSet(Vector_G->createDataSet(name_vector
            , PredType::NATIVE_DOUBLE, fspace_V));

        /*
         * Conversion Eigen -> DOUBLE
         */
        double rel_vel[dim_M];
        for (int i=0; i<dim_M; ++i){
            rel_vel[i] = lcp.q(i);
        }

        dataset_V->write(rel_vel, PredType::NATIVE_DOUBLE);

        /*
         * Close the dataset and the file.
         */        
        delete dataset_V;
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
}


}}} // namespace floe::lcp::solver


#endif // OPE_LCP_SOLVER_HPP