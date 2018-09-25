#include <iostream>
// #include <Eigen/SVD>
#include <limits>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/blas.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
// #include <boost/numeric/ublas/io.hpp>
#include <vector>
// using namespace Eigen;
#include <algorithm>

using namespace boost::numeric::ublas;

// #include <string>
// #include <H5Cpp.h>
// using namespace H5;

// const H5std_string file1( "matrix.h5" );

// using namespace std;

#include <string>
#include "H5Cpp.h"
using namespace H5;

// int main( int argc, char* argv[] )
// {

void blabla(int result);

int main (void)
{
	std::vector<int> v(10,0);
	v[3] = 10;
	std::cout << v[3] << "\n" << std::endl;

	std::cout << v.size() << "\n";

	int val[v.size()];
    for (size_t i=0; i<v.size(); ++i) {
    	val[i] = v[i];
    }

    std::cout << val[8] << std::endl;


    const H5std_string FILE_NAME( "io/outputs/selected_floes.h5" );
    const H5std_string GROUP_NAME_I( "selected_floe_ids" );
    H5File* file;
    file = new H5File( FILE_NAME, H5F_ACC_RDONLY );

    DataSet* dataset = new DataSet(file->openDataSet( GROUP_NAME_I ));
    /*
    * Get dataspace of the dataset.
    */
    DataSpace dataspace = dataset->getSpace();
    /*
    * Get the dimension size of each dimension in the dataspace and
    * display them.
    */
    hsize_t dims_out[1];
    dataspace.getSimpleExtentDims(dims_out);
    std::cout << "dim_out: " << dims_out[0] << "\n";
    int selection_tmp[dims_out[0]];
    dataset->read( selection_tmp, PredType::NATIVE_INT );

    std::cout << "selected_floe_ids: [\n";
    std::vector<std::size_t> selected_floe_ids;
    for (hsize_t i=0; i<dims_out[0]; ++i) {
        std::cout << selection_tmp[i] << ", ";
        selected_floe_ids.push_back(selection_tmp[i]);
    }
    std::cout << "]\n" << std::endl;

    delete dataset;
    delete file;
	// const H5std_string  FILE_NAME( "test" );
	// H5File* m_out_file;

	// /*
	//  * Try block to detect exceptions raised by any of the calls inside it
	//  */
	// try{
	//     /*
	//      * Turn off the auto-printing when failure occurs so that we can
	//      * handle the errors appropriately
	//      */
	//     Exception::dontPrint();
	//     /*
	//      * Create or Open a file.
	//      */
	//     try {
	//         m_out_file = new H5File( FILE_NAME, H5F_ACC_RDWR );
	//     } catch (...) {
	//         m_out_file = new H5File( FILE_NAME, H5F_ACC_TRUNC );

	//         /* write list of selected floes */
	//         hsize_t dim[2] = {1, v.size()};
	//         DataSpace space( 2, dim );

	//         int val[v.size()];
	//         for (size_t i=0; i<v.size(); ++i) {
	//         	val[i] = v[i];
	//         }

	//         std::cout << val << std::endl;
	//         DataSet* data_floes = new DataSet(m_out_file->createDataSet("selected_floe_ids", 
	//             PredType::NATIVE_INT, space));

	//         data_floes->write(val, PredType::NATIVE_INT);

	//         delete data_floes;
	//     }

	//     delete m_out_file;
	// }

	// // catch failure caused by the H5File operations
	// catch( FileIException error )
	// {
	// 	error.printError();
	// }
	// // catch failure caused by the DataSet operations
	// catch( DataSetIException error )
	// {
	// 	error.printError();
	// }

	// int vec_init[2] = {10,0};
	// int vec[4] = { vec_init[0], vec_init[1], 12, 12 };

	// std::cout << std::max(10,2);

	// int i=0;
	// std::cout << "vec: [";
	// for (i=0; i<4; i++){
	// 	std::cout << vec[i] << ", " << i << ", ";
	// 	if (i==2){break;}
	// } 
	// std::cout << ", " << i <<"]\n";
	// std::vector<int> v={6,0};

	// std::cout << v;
	// // std::vector<int>::iterator it=v.begin();
	/*for (std::size_t i=0; i<v.size(); ++i) {
		if (i==2 || i==4){ v[i] = 1;}
		else {v[i]=3;}
	}*/
 
    // std::vector<int>::iterator result = std::min_element(std::begin(v), std::end(v));
    // std::cout << "min element:" << *result << ", at: " << (result-v.begin()) << "\n";

    // blabla(*result);

	// int ll = -1;
	// int dd = 10;
	// int count = 0;
	// int bb = -1;

	// matrix< double > v(dd,dd);

	// std::vector<int> v(dd,0);

	// std::size_t dim=v.size();

	// for (int i=0; i<dim; ++i) {
 //        if (i == 3) {
 //            v[i] = 2;
 //        }
 //    }

	// std::vector<int>::iterator it;
	// it = std::find(v.begin(), v.end(), 2);
	// std::cout << "v(4) = " << v[4] << "\n";
	// std::cout << "element: " << *it << " find in: " << (it-v.begin()) << " find in related to the end: " << (it-v.end()) << "\n";

	// int index = (it-v.begin());
	// if (index==0) {std::cout << "sdfg\n";}

	// int i;
	// while (ll!=0 && dd==1 && bb!=0 && i<10){
	// 	int bb = 0;
	// 	std::cout << "je suis dans While depuis:" << count << "tour\n";

	// 	++i;

	// 	if (count>5){
	// 		dd = 1;
	// 		bb = 1;
	// 		++count;
	// 		continue;
	// 	}
	// 	if (count>11){
	// 		bb = 6;
	// 		ll = 0;
	// 	}
	// 	if (count>7){
	// 		bb = 0;
	// 	}
	// 	++count;
	// }

	// std::cout << "je suis sortie de While apres:" << count << "tour\n";
	// std::cout << "Mais non, je suis sortie de While apres:" << i << "tour\n";

	// int power = std::pow( 2, 20 );
	// int itermax = std::min(power , 10000);
	// std::cout << itermax << "\n";

	// double lcp_e = 1e-8;

	// bool solved=1;
	// const H5std_string FILE_NAME( "Select.h5" );
 //    // const H5std_string GROUP_NAME1( "Delassus Matrix of unsolved LCP" );
 //    // const H5std_string GROUP_NAME2( "Corresponding relative velocities" );
	// const int   MSPACE1_RANK = 1;   // Rank of the first dataset in memory
	// const int   MSPACE1_DIM = 50;   // Dataset size in memory
	// const int   MSPACE2_RANK = 1;   // Rank of the second dataset in memory
	// const int   MSPACE2_DIM = 4;    // Dataset size in memory
	// const int   FSPACE_RANK = 2;    // Dataset rank as it is stored in the file
	// const int   FSPACE_DIM1 = 8;    // Dimension sizes of the dataset as it is
	// const int   FSPACE_DIM2 = 12;   //  stored in the file
	// const int   MSPACE_RANK = 2;    // Rank of the first dataset in memory
	// const int   MSPACE_DIM1 = 8;    // We will read dataset back from the file
	// const int   MSPACE_DIM2 = 9;    //  to the dataset in memory with these
	//                 //  dataspace parameters
	// const int   NPOINTS = 4;    // Number of points that will be selected
	//                 //  and overwritten

	// bool a = (0==0)? 1:0;

	// const int test_idx = 2;
	// const int solver_used = 3;

	// int size_delassus = 12;
	// MatrixXd BMB(size_delassus,size_delassus);
	// BMB << 0.00152306,0.00167011,0.00159376,-0.00026744,0.00026744,-0.000166789,0.000166789,-0.000162222,0.000162222,0,0,0,
	// 0.00167011,0.00241552,0.00202767,-0.00113476,0.00113476,-0.00102235,0.00102235,-0.000999151,0.000999151,0,0,0,
	// 0.00159376,0.00202767,0.00180468,-0.000754591,0.000754591,-0.000648139,0.000648139,-0.0006348,0.0006348,0,0,0,
	// -0.00026744,-0.00113476,-0.000764591,0.00314846,-0.00314846,0.00312383,-0.00312383,0.00310109,-0.00310109,1,0,0,
	// 0.00026744,0.00113476,0.000754591,-0.00314846,0.00314846,-0.00312383,0.00312383,-0.00310109,0.00310109,1,0,0,
	// -0.000166789,-0.00102235,-0.000648139,0.00312383,-0.00312383,0.00310586,-0.00310586,0.00308348,-0.00308348,0,1,0,
	// 0.000166789,0.00102235,0.000648139,-0.00312383,0.00312383,-0.00310586,0.00310586,-0.00308348,0.00308348,0,1,0,
	// -0.000162222,-0.000999151,-0.0006348,0.00310109,-0.00310109,0.00308348,-0.00308348,0.00306168,-0.00306168,0,0,1,
	// 0.000162222,0.000999151,0.0006348,-0.00310109,0.00310109,-0.00308348,0.00308348,-0.00306168,0.00306168,0,0,1,
	// 0.7,0,0,-1,-1,0,0,0,0,0,0,0,
	// 0,0.7,0,0,0,-1,-1,0,0,0,0,0,
	// 0,0,0.7,0,0,0,0,-1,-1,0,0,0;

	// BMB(0,2) = 12;

	// MatrixXd BB(9,9);
	// BB << 0.00152306,0.00167011,0.00159376,-0.00026744,0.00026744,-0.000166789,0.000166789,-0.000162222,0.000162222,
	// 0.00167011,0.00241552,0.00202767,-0.00113476,0.00113476,-0.00102235,0.00102235,-0.000999151,0.000999151,
	// 0.00159376,0.00202767,0.00180468,-0.000754591,0.000754591,-0.000648139,0.000648139,-0.0006348,0.0006348,
	// -0.00026744,-0.00113476,-0.000754591,0.00314846,-0.00314846,0.00312383,-0.00312383,0.00310109,-0.00310109,
	// 0.00026744,0.00113476,0.000754591,-0.00314846,0.00314846,-0.00312383,0.00312383,-0.00310109,0.00310109,
	// -0.000166789,-0.00102235,-0.000648139,0.00312383,-0.00312383,0.00310586,-0.00310586,0.00308348,-0.00308348,
	// 0.000166789,0.00102235,0.000648139,-0.00312383,0.00312383,-0.00310586,0.00310586,-0.00308348,0.00308348,
	// -0.000162222,-0.000999151,-0.0006348,0.00310109,-0.00310109,0.00308348,-0.00308348,0.00306168,-0.00306168,
	// 0.000162222,0.000999151,0.0006348,-0.00310109,0.00310109,-0.00308348,0.00308348,-0.00306168,0.00306168;

	// // const H5std_string FILE_NAME("/Users/matthiasrabatel/Travail/outputs_mycode/matrix.h5");
 //    const H5std_string GROUP_NAME_I( "solved" ); // root group
 //    const H5std_string GROUP_NAME_II( "unsolved" ); // root group
 //    const H5std_string GROUP_NAME1( "M" );
 //    const H5std_string GROUP_NAME2( "q" );
 //    const H5std_string GROUP_NAME3( "z" );
 //    const H5std_string LCP_error( "LCP error" );
 //    const H5std_string Last_Memb( "Last LCP" );
 //    const H5std_string Idx_solver( "Which solver" ); // Information on which solver and 
 //    H5std_string GROUP_TEMP;

 //    if (solved){
 //        GROUP_TEMP = "solved";
 //    }
 //    else {
 //        GROUP_TEMP = "unsolved";
 //    }

 //    /*
 //     * Try block to detect exceptions raised by any of the calls inside it
 //     */
 //    // try{
 //        /*
 //         * Turn off the auto-printing when failure occurs so that we can
 //         * handle the errors appropriately
 //         */
 //        // Exception::dontPrint();
 //        /*
 //         * Create or Open a file.
 //         */
 //        H5File* file;
 //        try {
 //            file = new H5File( FILE_NAME, H5F_ACC_RDWR );
 //        } catch (...) {
 //            file = new H5File( FILE_NAME, H5F_ACC_TRUNC );
 //            Group* M_solved = new Group(file->createGroup(GROUP_NAME_I));
 //            Group(M_solved->createGroup(GROUP_NAME1));
 //            Group(M_solved->createGroup(GROUP_NAME2));
 //            Group(M_solved->createGroup(GROUP_NAME3));

 //            Group* M_unsolved = new Group(file->createGroup(GROUP_NAME_II));
 //            Group(M_unsolved->createGroup(GROUP_NAME1));
 //            Group(M_unsolved->createGroup(GROUP_NAME2));
 //            Group(M_unsolved->createGroup(GROUP_NAME3));

 //            hsize_t dim_LM[1] = {1};
 //            DataSpace space_LM( 1, dim_LM );
 //            DataSet(M_solved->createDataSet( Last_Memb, PredType::NATIVE_INT, space_LM ));
 //            DataSet(M_unsolved->createDataSet( Last_Memb, PredType::NATIVE_INT, space_LM ));
            
 //            hsize_t dim_idx_solver[2] = {1, 2};
 //            hsize_t maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED}; // unlimited dataspace
 //            DataSpace space_solver( 2, dim_idx_solver, maxdims );
            
 //            DSetCreatPropList prop; // Modify dataset creation property to enable chunking
 //            hsize_t chunk_dims[2] = {1, 2}; // with extendible dataset we cannot use contiguous but chunked dataset
 //            prop.setChunk(2, chunk_dims);
 //            std::cout << "je suis bien arrive ici: -1\n";
 //            DataSet(M_solved->createDataSet( Idx_solver, PredType::NATIVE_INT, space_solver, prop ));
 //            DataSet(M_unsolved->createDataSet( Idx_solver, PredType::NATIVE_INT, space_solver, prop ));

 //            hsize_t maxdims_le[1] = {H5S_UNLIMITED};
 //            DataSpace space_LE( 1, dim_LM, maxdims_le );
 //            DSetCreatPropList prop_le; // Modify dataset creation property to enable chunking
 //            hsize_t chunk_dims_le[1] = {1}; // with extendible dataset we cannot use contiguous but chunked dataset
 //            prop_le.setChunk(1, chunk_dims_le);
 //            DataSet(M_solved->createDataSet( LCP_error, PredType::NATIVE_DOUBLE, space_LE, prop_le ));
 //            DataSet(M_unsolved->createDataSet( LCP_error, PredType::NATIVE_DOUBLE, space_LE, prop_le ));

	// 		std::cout << "je suis bien arrive ici: 0\n";

 //            delete M_unsolved;
 //            delete M_solved;
 //        }
	// 	/*
 //         * Recovering the number of LCP
 //         */
 //        Group *Matrix_G, *Vector_G, *Z_G, *Root;

 //        Group* M_solved = new Group(file->openGroup(GROUP_NAME_I));
 //        Group* MS = new Group(M_solved->openGroup(GROUP_NAME1));
 //        hsize_t nb_lcp_sol = MS->getNumObjs();

 //        Group* M_unsolved = new Group(file->openGroup(GROUP_NAME_II));
 //        Group* MU = new Group(M_unsolved->openGroup(GROUP_NAME1));
 //        hsize_t nb_lcp_unsol = MU->getNumObjs();

 //        hsize_t nb_lcp = nb_lcp_sol + nb_lcp_unsol;

 //        delete M_unsolved;
 //        delete M_solved;
 //        delete MS;
 //        delete MU;

 //        Root = new Group(file->openGroup(GROUP_TEMP));
 //        Matrix_G = new Group(Root->openGroup(GROUP_NAME1));
 //        Vector_G = new Group(Root->openGroup(GROUP_NAME2));
 //        Z_G = new Group(Root->openGroup(GROUP_NAME3));

 //        hsize_t nb_lcp_temp = Matrix_G->getNumObjs();

 //        bool G_exist = 1;
 //        if (nb_lcp_temp==0) {
 //        	G_exist=0;
 //        	DataSet* dataset_solver = new DataSet(Root->openDataSet( Idx_solver ));
 //        	int idx_solv[2] = { test_idx , solver_used };
 //        	std::cout << "idx_solv: " << idx_solv[0] << idx_solv[1] << "\n";
 //        	dataset_solver->write(idx_solv, PredType::NATIVE_INT);


 //        	DataSet* dataset_LE = new DataSet(Root->openDataSet( LCP_error ));
 //            // double lcp_e = lcp_err;
 //            dataset_LE->write(&lcp_e, PredType::NATIVE_DOUBLE);
 //            std::cout << "coucou: 1\n";

 //        }

 //        /*
 //         * Comparison to the previous LCP failure (to prevent similar LCP)
 //         */
 //        bool isnt_same_LCP = 1;

 //        if (G_exist) {
 //            int last_lcp[1];
 //            DataSet* dataset_LM = new DataSet(Root->openDataSet( Last_Memb ));
 //            dataset_LM->read( last_lcp, PredType::NATIVE_INT );


 //            const H5std_string name_data_pre = std::to_string(last_lcp[0]);

 //            DataSet* dataset_pre = new DataSet(Matrix_G->openDataSet( name_data_pre ));

 //            DataSpace fspace1 = dataset_pre->getSpace();
 //            std::size_t dim_out = std::sqrt( fspace1.getSelectNpoints() );

   
 //            if (dim_out==size_delassus){
 //                double data_out[dim_out][dim_out];

 //                dataset_pre->read( data_out, PredType::NATIVE_DOUBLE ); 

	//          	/*
	//          	 * Check if matrix already exists? (A large number of attempt to solve LCP)
	//          	 */
	//         	MatrixXd Diff( dim_out , dim_out );
	//         	for (int i=0; i<dim_out; ++i){
	//         		for (int j=0; j<dim_out; ++j){
	//         			const double val_rel = std::min( std::abs(BMB(i,j)) , std::abs(data_out[i][j]) );
	//         			double div = val_rel;
	//         			if (val_rel==0) {div = 1.0;}
	//         			const double val_rel_a = (BMB(i,j) - data_out[i][j])/div;

	//         			Diff(i,j) = std::max( std::abs( val_rel_a ) , 0.0);
	//         		}
	//         	}
	//         	isnt_same_LCP = Diff.norm() > 1e-7;

	//         	std::cout << "it's ok, is_same:" << isnt_same_LCP << "\n";
 //            } 
 //        }

 //        /*
 //         * Create dataspace for the dataset in the file.
 //         */
 //        const H5std_string name_matrix = std::to_string(nb_lcp+1);

 //        if (isnt_same_LCP) {
 //        	std::cout << "loop is_same_LCP=0 \n";
	//         hsize_t dim_space_M[2];
	//         dim_space_M[0] = size_delassus;
	//         dim_space_M[1] = size_delassus;
	//         DataSpace fspace_M( 2, dim_space_M );
	//         /*
	//          * Create dataset and write it into the file.
	//          */
	//         DataSet* dataset_M = new DataSet(Matrix_G->createDataSet(name_matrix
	//             , PredType::NATIVE_DOUBLE, fspace_M));

	// 		/*
	//     	 * Conversion Eigen -> DOUBLE
	//     	 */
	//     	double Delassus[size_delassus][size_delassus];
	//     	for (int i=0; i<size_delassus; ++i){
	//     		for (int j=0; j<size_delassus; ++j){
	// 				Delassus[i][j] = BMB(i,j);
	//     		}
	//     	}

	//         dataset_M->write(Delassus, PredType::NATIVE_DOUBLE);

	//         /*
	//          * Close the dataset and the file.
	//          */        
	//         delete dataset_M;

	//         DataSet* dataset_LM = new DataSet(Root->openDataSet( Last_Memb ));
 //            const int nb_LM[1] = {static_cast<int>(nb_lcp+1)};
 //            dataset_LM->write(nb_LM, PredType::NATIVE_INT);
 //            delete dataset_LM;

 //            /*
 //             * Save information on solvers 
 //             */
 //            if (G_exist){
 //            	std::cout << "je suis bien arrive ici: 1\n";
	//             DataSet* dataset_solver = new DataSet(Root->openDataSet( Idx_solver ));
	//             DataSpace space_solver = dataset_solver->getSpace();
	//             hsize_t dim_curr[2]; // dimension of the dataset
	//             space_solver.getSimpleExtentDims( dim_curr, NULL); // retrieves the current dimensions 
	//             std::cout << "dim_curr: " << dim_curr[0] << dim_curr[1] << "\n";
	//             hsize_t ext_size[2] = { dim_curr[0]+1, dim_curr[1]}; 
	//             std::cout << "ext_size: " << ext_size[0] << ext_size[1] << "\n";
	//             dataset_solver->extend( ext_size ); // extension with one new line 
	  
	//             std::cout << "je suis bien arrive ici: 2\n";
	//             DataSpace fspace2 = dataset_solver->getSpace();
	//             hsize_t dim2[2] = {1,2}; 
	//             hsize_t offset2[2] = {dim_curr[0], 0};
	//             fspace2.selectHyperslab( H5S_SELECT_SET, dim2, offset2); // selection of the hyperslab
	//             DataSpace mspace2( 2, dim2 );
	//             int idx_solv[2] = { test_idx , solver_used };
	//             std::cout << "je suis bien arrive ici: 3\n";
	//             std::cout << "idx_solv: " << idx_solv[0] << idx_solv[1] << "\n";
	//             std::cout << "offset2: " << offset2[0] << offset2[1] << "\n";
	//             std::cout << "dim2: " << dim2[0] << dim2[1] << "\n";
	//             dataset_solver->write(idx_solv, PredType::NATIVE_INT, mspace2, fspace2); // write in the hyperslab
	//             delete dataset_solver;
	//             std::cout << "je suis bien arrive ici: 4\n";


	//             /*
 //                 * Save information on LCP error with extendible dataset
 //                 */                
 //                DataSet* dataset_LE = new DataSet(Root->openDataSet( LCP_error ));
 //                DataSpace space_LE = dataset_LE->getSpace();
 //                hsize_t dim_curr_le[1]; // dimension of the dataset
 //                space_LE.getSimpleExtentDims( dim_curr_le, NULL); // retrieves the current dimensions 
 //                hsize_t ext_size_le[1] = { dim_curr_le[0]+1}; 
 //                dataset_LE->extend( ext_size_le ); // extension with one new line 
      
 //                DataSpace fspace_le = dataset_LE->getSpace();
 //                hsize_t dim_le[1] = {1}; 
 //                hsize_t offset_le[1] = {dim_curr_le[0]};
 //                fspace_le.selectHyperslab( H5S_SELECT_SET, dim_le, offset_le); // selection of the hyperslab
 //                DataSpace mspace_le( 1, dim_le );
 //                // double lcp_e = lcp_err;
 //                std::cout << "dim_le: " << dim_le[0] << "\n";
 //                std::cout << "je suis bien arrive ici: 5\n";
 //                dataset_LE->write(&lcp_e, PredType::NATIVE_DOUBLE, mspace_le, fspace_le); // write in the hyperslab

 //                delete dataset_LE;
 //            }
            
	//         /*-----------------------------------------------------------------------------------------
	//          * new dataset for relative velocities
	//          *---------------------------------------------------------------------------------------*/
	//         const H5std_string name_vector = std::to_string(nb_lcp+1);
	//         /*
	//          * Create dataspace for the dataset in the file.
	//          */

	//         hsize_t dim_space_V[1];
	//         dim_space_V[0] = size_delassus;

	//         DataSpace fspace_V( 1, dim_space_V );
	//         /*
	//          * Create dataset and write it into the file.
	//          */

	//         DataSet* dataset_V = new DataSet(Vector_G->createDataSet(name_vector
	//             , PredType::NATIVE_DOUBLE, fspace_V));

	//         /*
	//          * Conversion Eigen -> DOUBLE
	//          */
	//         double rel_vel[size_delassus];
	//         for (int i=0; i<size_delassus; ++i){
	//             rel_vel[i] = BMB(i,1);
	//         }

	//         dataset_V->write(rel_vel, PredType::NATIVE_DOUBLE);

	//        /*
	//         * Close the dataset and the file.
	//         */        
	//         delete dataset_V;

	//         const H5std_string name_z = name_matrix;
 //            /*
 //             * Create dataspace for the dataset in the file.
 //             */

 //            hsize_t dim_z[1];
 //            dim_z[0] = size_delassus;
 //            DataSpace fspace_z( 1, dim_z );
 //            /*
 //             * Create dataset and write it into the file.
 //             */
 //            DataSet* dataset_z = new DataSet(Z_G->createDataSet(name_z
 //                , PredType::NATIVE_DOUBLE, fspace_z));

 //            /*
 //             * Conversion Eigen -> DOUBLE
 //             */
 //            double lcp_z[size_delassus];
 //            for (int i=0; i<size_delassus; ++i){
 //                lcp_z[i] = BMB(i,1);
 //            }
 //            dataset_z->write(lcp_z, PredType::NATIVE_DOUBLE);

 //            /*
 //             * Close the dataset and the file.
 //             */        
 //            delete dataset_z;

 //    	}

 //        delete file;

        return 0;
}

void blabla(int result) {
	std::cout << "I'm able to return: " << result << "\n";
}

    // /*
    //  * Try block to detect exceptions raised by any of the calls inside it
    //  */
    // try{
	   //  /*
	   //   * Turn off the auto-printing when failure occurs so that we can
	   //   * handle the errors appropriately
	   //   */
	   //  Exception::dontPrint();
	    
	   //   * Create or Open a file.
	     
	   //  H5File* file;
	   //  if (!a){
	   //  	file = new H5File( FILE_NAME, H5F_ACC_TRUNC );	
	   //  	std::cout << "new!!\n";
	   //  } 
	   //  else {
	   //  	file = new H5File( FILE_NAME, H5F_ACC_RDWR );
	   //  	std::cout << "deja cree\n";
	   //  }
   // }  // end of try block
   // // catch failure caused by the H5File operations
   // catch( FileIException error )
   // {
   //  error.printErrorStack();
   // }
   // // catch failure caused by the DataSet operations
   // catch( DataSetIException error )
   // {
   //  error.printErrorStack();
   // }
   // // catch failure caused by the DataSpace operations
   // catch( DataSpaceIException error )
   // {
   //  error.printErrorStack();
   // }


/*/////////////////////////////////////////////////////////////////////////////////////////////////*/
		/*
	    * Create property list for a dataset and set up fill values.
	    */
	 //    int fillvalue = 0;   /* Fill value for the dataset */
	 //    DSetCreatPropList plist;
	 //    plist.setFillValue(PredType::NATIVE_INT, &fillvalue);

	 //    hsize_t fdim[] = {FSPACE_DIM1, FSPACE_DIM2}; // dim sizes of ds (on disk)
	 //    DataSpace fspace( FSPACE_RANK, fdim );

	 //    DataSet* dataset = new DataSet(file->createDataSet(
	 //        DATASET_NAME, PredType::NATIVE_INT, fspace, plist));

		// hsize_t start[2]; // Start of hyperslab
	 //    hsize_t stride[2]; // Stride of hyperslab
	 //    hsize_t count[2];  // Block count
	 //    hsize_t block[2];  // Block sizes
	 //    start[0]  = 0; start[1]  = 1;
	 //    stride[0] = 4; stride[1] = 3;
	 //    count[0]  = 2; count[1]  = 4;
	 //    block[0]  = 3; block[1]  = 2;
	 //    fspace.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);
	 //    /*
	 //     * Create dataspace for the first dataset.
	 //     */
	 //    hsize_t dim1[] = {MSPACE1_DIM};  /* Dimension size of the first dataset
	 //                                       (in memory) */
	 //    DataSpace mspace1( MSPACE1_RANK, dim1 );
	 //    /*
	 //     * Select hyperslab.
	 //     * We will use 48 elements of the vector buffer starting at the
	 //     * second element.  Selected elements are 1 2 3 . . . 48
	 //     */
	 //    start[0]  = 1;
	 //    stride[0] = 1;
	 //    count[0]  = 48;
	 //    block[0]  = 1;
	 //    mspace1.selectHyperslab( H5S_SELECT_SET, count, start, stride, block);

	 //    int    vector[MSPACE1_DIM]; // vector buffer for dset
	 //    /*
	 //     * Buffer initialization.
	 //     */
	 //    vector[0] = vector[MSPACE1_DIM - 1] = -1;
	 //    for (i = 1; i < MSPACE1_DIM - 1; i++)
	 //        vector[i] = i;
	 //    dataset->write( vector, PredType::NATIVE_INT, mspace1, fspace );
	 //    /*
	 //     * Reset the selection for the file dataspace fid.
	 //     */
	 //    fspace.selectNone();
	 //    /*
	 //     * Create dataspace for the second dataset.
	 //     */
	 //    hsize_t dim2[] = {MSPACE2_DIM};  /* Dimension size of the second dataset
	 //                                       (in memory */
	 //    DataSpace mspace2( MSPACE2_RANK, dim2 );
	 //    /*
	 //     * Select sequence of NPOINTS points in the file dataspace.
	 //     */
	 //    hsize_t coord[NPOINTS][FSPACE_RANK]; /* Array to store selected points
	 //                                            from the file dataspace */
	 //    coord[0][0] = 0; coord[0][1] = 0;
	 //    coord[1][0] = 3; coord[1][1] = 3;
	 //    coord[2][0] = 3; coord[2][1] = 5;
	 //    coord[3][0] = 5; coord[3][1] = 6;
	 //    fspace.selectElements( H5S_SELECT_SET, NPOINTS, (const hsize_t *)coord);
	 //    /*
	 //     * Write new selection of points to the dataset.
	 //     */
	 //    int    values[] = {53, 59, 61, 67};  /* New values to be written */
	 //    dataset->write( values, PredType::NATIVE_INT, mspace2, fspace );
	     
	 //    /*
	 //     * Close the dataset and the file.
	 //     */
	 //    delete dataset;
	 //    delete file;

  //   	/*
	 //     * Open the file.
	 //     */

	 //    file = new H5File( FILE_NAME, H5F_ACC_RDONLY );
	 //    /*
	 //     * Open the dataset.
	 //     */
	 //    dataset = new DataSet( file->openDataSet( DATASET_NAME ));
	 //    /*
	 //     * Get dataspace of the dataset.
	 //     */
	 //    fspace = dataset->getSpace();

	 //    start[0] = 1; start[1] = 2;
	 //    block[0] = 1; block[1] = 1;
	 //    stride[0] = 1; stride[1] = 1;
	 //    count[0]  = 3; count[1]  = 4;
	 //    fspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);

	 //    start[0] = 2; start[1] = 4;
	 //    block[0] = 1; block[1] = 1;
	 //    stride[0] = 1; stride[1] = 1;
	 //    count[0]  = 6; count[1]  = 5;
	 //    fspace.selectHyperslab(H5S_SELECT_OR, count, start, stride, block);
	 //    /*
	 //     * Create memory dataspace.
	 //     */
	 //    hsize_t mdim[] = {MSPACE_DIM1, MSPACE_DIM2}; /* Dimension sizes of the
	 //                                                   dataset in memory when we
	 //                                                   read selection from the
	 //                                                   dataset on the disk */
	 //    DataSpace mspace(MSPACE_RANK, mdim);
	 //    /*
	 //     * Select two hyperslabs in memory. Hyperslabs has the same
	 //     * size and shape as the selected hyperslabs for the file dataspace.
	 //     */
	 //    start[0] = 0; start[1] = 0;
	 //    block[0] = 1; block[1] = 1;
	 //    stride[0] = 1; stride[1] = 1;
	 //    count[0]  = 3; count[1]  = 4;
	 //    mspace.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
	 //    start[0] = 1; start[1] = 2;
	 //    block[0] = 1; block[1] = 1;
	 //    stride[0] = 1; stride[1] = 1;
	 //    count[0]  = 6; count[1]  = 5;
	 //    mspace.selectHyperslab(H5S_SELECT_OR, count, start, stride, block);
	 //    /*
	 //     * Initialize data buffer.
	 //     */
	 //    int matrix_out[MSPACE_DIM1][MSPACE_DIM2];
	 //    for (i = 0; i < MSPACE_DIM1; i++)
	 //        for (j = 0; j < MSPACE_DIM2; j++)
	 //        matrix_out[i][j] = 0;
	 //    /*
	 //     * Read data back to the buffer matrix.
	 //     */
	 //    dataset->read(matrix_out, PredType::NATIVE_INT, mspace, fspace);

	 //    for (i=0; i < MSPACE_DIM1; i++)
	 //    {
	 //        for(j=0; j < MSPACE_DIM2; j++)
	 //        cout << matrix_out[i][j] << "  ";
	 //        cout << endl;
	 //    }
	 //    /*
	 //     * Close the dataset and the file.
	 //     */
	 //    delete dataset;
	 //    delete file;

/*/////////////////////////////////////////////////////////////////////////////////////////////////*/

// 	    Group* Matrix_G;
// 		try {
//         	Matrix_G = new Group(file->openGroup(GROUP_NAME));
//     	} catch (...) {
// 	        /* Create group for floe shapes */
// 	       	Matrix_G = new Group(file->createGroup(GROUP_NAME));
// 	    }

// 	    hsize_t Group_size = Matrix_G->getNumObjs();
// 		const H5std_string name_matrix = std::to_string(Group_size+1);
// 	    /*
// 	     * Create dataspace for the dataset in the file.
// 	     */
// 	    const int dim_DM1 = BMB.rows();
// 	    const int dim_DM2 = BMB.cols();

// 	    hsize_t dim_space[2] = {dim_DM1+9, dim_DM2+9}; // dim sizes of ds (on disk)
// 	    DataSpace fspace( 2, dim_space );
// 	    /*
// 	     * Create dataset and write it into the file.
// 	     */

//     	DataSet* dataset = new DataSet(Matrix_G->createDataSet(name_matrix
// 			, PredType::NATIVE_DOUBLE, fspace));

//     	/*
//     	 * Conversion Eigen -> DOUBLE
//     	 */
//     	double Delassus[dim_DM1][dim_DM2];
//     	for (int i=0; i<dim_DM1; ++i){
//     		for (int j=0; j<dim_DM2; ++j){
// 				Delassus[i][j] = BMB(i,j);
//     		}
//     	}

//     // 	for (int i=0; i<9; ++i){
//     // 		for (int j=0; j<9; ++j){
// 				// Delassus[i+12][j+12] = BB(i,j);
//     // 		}
//     // 	}

//     	dataset->write(Delassus, PredType::NATIVE_DOUBLE);

// 	    /*
// 	     * Close the dataset and the file.
// 	     */
// 	    delete dataset;
// 	    delete file;	    
//    	}  // end of try block
//    	// catch failure caused by the H5File operations

// 	catch( FileIException error ){
// 	    error.printErrorStack();
// 	    return -1;
//    	}
//    	// catch failure caused by the DataSet operations
//    	catch( DataSetIException error )
//    	{
//     	error.printErrorStack();
//     	return -1;
//    	}
//    	// catch failure caused by the DataSpace operations
//    	catch( DataSpaceIException error )
//    	{
//     	error.printErrorStack();
//     	return -1;
//    	}
//    	return 0;
// }

	// double blabla[2] = {0.0,0.0};
	// for (int i=0; i<10; ++i){
	// 	blabla[1] += 2;
	// }

	// double qsdf = (blabla[0]==0)? 10 : 100/blabla[0];
	// std::cout << qsdf <<"\n";
	

    // int size_delassus = 4;
    // Eigen::MatrixXd BMB(size_delassus,size_delassus);
    // BMB(0,0) = 7, BMB(0,1) = 1, BMB(0,2) = 11, BMB(0,3) = 10;
    // BMB(1,0) = 2, BMB(1,1) = 6, BMB(1,2) = 5, BMB(1,3) = 2;
    // BMB(2,0) = 8, BMB(2,1) = 11, BMB(2,2) = 3, BMB(2,3) = 8;
    // BMB(3,0) = 6, BMB(3,1) = 9, BMB(3,2) = 3, BMB(3,3) = 6;

	// double lcp_failed_stats[1] = {0.0};

	// for (int i=0;i<1e2; ++i){

	// 	// MatrixXd AA(12,12);
	// 	// AA << 0.00152306,0.00167011,0.00159376,-0.00026744,0.00026744,-0.000166789,0.000166789,-0.000162222,0.000162222,0,0,0,
	// 	// 0.00167011,0.00241552,0.00202767,-0.00113476,0.00113476,-0.00102235,0.00102235,-0.000999151,0.000999151,0,0,0,
	// 	// 0.00159376,0.00202767,0.00180468,-0.000754591,0.000754591,-0.000648139,0.000648139,-0.0006348,0.0006348,0,0,0,
	// 	// -0.00026744,-0.00113476,-0.000754591,0.00314846,-0.00314846,0.00312383,-0.00312383,0.00310109,-0.00310109,1,0,0,
	// 	// 0.00026744,0.00113476,0.000754591,-0.00314846,0.00314846,-0.00312383,0.00312383,-0.00310109,0.00310109,1,0,0,
	// 	// -0.000166789,-0.00102235,-0.000648139,0.00312383,-0.00312383,0.00310586,-0.00310586,0.00308348,-0.00308348,0,1,0,
	// 	// 0.000166789,0.00102235,0.000648139,-0.00312383,0.00312383,-0.00310586,0.00310586,-0.00308348,0.00308348,0,1,0,
	// 	// -0.000162222,-0.000999151,-0.0006348,0.00310109,-0.00310109,0.00308348,-0.00308348,0.00306168,-0.00306168,0,0,1,
	// 	// 0.000162222,0.000999151,0.0006348,-0.00310109,0.00310109,-0.00308348,0.00308348,-0.00306168,0.00306168,0,0,1,
	// 	// 0.7,0,0,-1,-1,0,0,0,0,0,0,0,
	// 	// 0,0.7,0,0,0,-1,-1,0,0,0,0,0,
	// 	// 0,0,0.7,0,0,0,0,-1,-1,0,0,0;

	// 	// MatrixXd BB(9,9);
	// 	// BB << 0.00152306,0.00167011,0.00159376,-0.00026744,0.00026744,-0.000166789,0.000166789,-0.000162222,0.000162222,
	// 	// 0.00167011,0.00241552,0.00202767,-0.00113476,0.00113476,-0.00102235,0.00102235,-0.000999151,0.000999151,
	// 	// 0.00159376,0.00202767,0.00180468,-0.000754591,0.000754591,-0.000648139,0.000648139,-0.0006348,0.0006348,
	// 	// -0.00026744,-0.00113476,-0.000754591,0.00314846,-0.00314846,0.00312383,-0.00312383,0.00310109,-0.00310109,
	// 	// 0.00026744,0.00113476,0.000754591,-0.00314846,0.00314846,-0.00312383,0.00312383,-0.00310109,0.00310109,
	// 	// -0.000166789,-0.00102235,-0.000648139,0.00312383,-0.00312383,0.00310586,-0.00310586,0.00308348,-0.00308348,
	// 	// 0.000166789,0.00102235,0.000648139,-0.00312383,0.00312383,-0.00310586,0.00310586,-0.00308348,0.00308348,
	// 	// -0.000162222,-0.000999151,-0.0006348,0.00310109,-0.00310109,0.00308348,-0.00308348,0.00306168,-0.00306168,
	// 	// 0.000162222,0.000999151,0.0006348,-0.00310109,0.00310109,-0.00308348,0.00308348,-0.00306168,0.00306168;

	// 	int size_delassus = 3*AA.rows()/4;
 //        Eigen::MatrixXd BMB(size_delassus,size_delassus);
 //        for (int i=0;i<size_delassus;++i){
 //            for (int j=0;j<size_delassus;++j){
 //                BMB(i,j) = AA(i,j);
 //            }
 //        }

	// 	JacobiSVD<MatrixXd> svd(BMB);
	//     double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
	    
	//     lcp_failed_stats[0] += cond;
	// }

 //    std::cout << "blabla:" << lcp_failed_stats[0] << "\n";

    // Optimisation:
 //    MatrixXd Id = MatrixXd::Identity(12,12);
 //    long coef = 1e2;
 //    MatrixXd AA_tilde(12,12);
 //    AA_tilde = coef*Id*AA*Id*coef;

	// JacobiSVD<MatrixXd> svd2(AA_tilde);
 //    cond = svd2.singularValues()(0) / svd2.singularValues()(svd2.singularValues().size()-1);
    
 //    std::cout << "After optimisation, the condition number is:" << cond << "\n";
 //    std::cout << "BB_tilde is:" << AA_tilde << "\n";

    // Eigen::JacobiSVD<Eigen::MatrixXd> svd(AA);
    // double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
    
    // std::cout << "condition number is:" << cond << "\n";

    // long a(4);
    // double i(4.5969);

    // std::cout << "blabla:" << i/a << "\n";
    // std::cout << "Maximum value for int: " << std::numeric_limits<int>::max() << '\n';
    // std::cout << "Maximum value for double: " << std::numeric_limits<double>::max() << '\n';
    // std::cout << "Maximum value for long: " << std::numeric_limits<long>::max() << '\n';