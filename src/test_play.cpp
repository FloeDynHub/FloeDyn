#include <iostream>
// #include <Eigen/SVD>
#include <limits>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/blas.hpp>
#include <boost/multi_array.hpp>

// using namespace Eigen;

// #include <string>
// #include <H5Cpp.h>
// using namespace H5;

// const H5std_string file1( "matrix.h5" );

using std::cout;
using std::endl;

#include <string>
#include "H5Cpp.h"
using namespace H5;


// int main( int argc, char* argv[] )
// {


int main (void)
{

	int power = std::pow( 2, 20 );
	int itermax = std::min(power , 10000);
	std::cout << itermax << "\n";


 //    const H5std_string FILE_NAME( "Select.h5" );
 //    const H5std_string GROUP_NAME1( "Delassus Matrix of unsolved LCP" );
 //    const H5std_string GROUP_NAME2( "Corresponding relative velocities" );
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

	// int size_delassus = 12;
	// MatrixXd BMB(size_delassus,size_delassus);
	// BMB << 0.00152306,0.00167011,0.00159376,-0.00026744,0.00026744,-0.000166789,0.000166789,-0.000162222,0.000162222,0,0,0,
	// 0.00167011,0.00241552,0.00202767,-0.00113476,0.00113476,-0.00102235,0.00102235,-0.000999151,0.000999151,0,0,0,
	// 0.00159376,0.00202767,0.00180468,-0.000754591,0.000754591,-0.000648139,0.000648139,-0.0006348,0.0006348,0,0,0,
	// -0.00026744,-0.00113476,-0.000754591,0.00314846,-0.00314846,0.00312383,-0.00312383,0.00310109,-0.00310109,1,0,0,
	// 0.00026744,0.00113476,0.000754591,-0.00314846,0.00314846,-0.00312383,0.00312383,-0.00310109,0.00310109,1,0,0,
	// -0.000166789,-0.00102235,-0.000648139,0.00312383,-0.00312383,0.00310586,-0.00310586,0.00308348,-0.00308348,0,1,0,
	// 0.000166789,0.00102235,0.000648139,-0.00312383,0.00312383,-0.00310586,0.00310586,-0.00308348,0.00308348,0,1,0,
	// -0.000162222,-0.000999151,-0.0006348,0.00310109,-0.00310109,0.00308348,-0.00308348,0.00306168,-0.00306168,0,0,1,
	// 0.000162222,0.000999151,0.0006348,-0.00310109,0.00310109,-0.00308348,0.00308348,-0.00306168,0.00306168,0,0,1,
	// 0.7,0,0,-1,-1,0,0,0,0,0,0,0,
	// 0,0.7,0,0,0,-1,-1,0,0,0,0,0,
	// 0,0,0.7,0,0,0,0,-1,-1,0,0,0;

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

	// try{
 //        /*
 //         * Turn off the auto-printing when failure occurs so that we can
 //         * handle the errors appropriately
 //         */
 //        Exception::dontPrint();
 //        /*
 //         * Create or Open a file.
 //         */
 //        H5File* file;
 //        try {
 //            file = new H5File( FILE_NAME, H5F_ACC_RDWR );
 //        } catch (...) {
 //            file = new H5File( FILE_NAME, H5F_ACC_TRUNC );
 //        }
        
 //        /*
 //         * Create or Open the groups
 //         */
 //        Group* Matrix_G;
 //        Group* Vector_G;
 //        bool G_exist;
 //        try {
 //            Matrix_G = new Group(file->openGroup(GROUP_NAME1));
 //            Vector_G = new Group(file->openGroup(GROUP_NAME2));
 //            G_exist = 1;
 //        } catch (...) {
 //            /* Create group for floe shapes */
 //            Matrix_G = new Group(file->createGroup(GROUP_NAME1));
 //            Vector_G = new Group(file->createGroup(GROUP_NAME2));
 //            G_exist = 0;
 //        }

 //        hsize_t Group_size = Matrix_G->getNumObjs();
 //        const H5std_string name_matrix = std::to_string(Group_size+1);

	//     const int dim_M = size_delassus;
 //        /*
 //         * Conversion Eigen -> DOUBLE
 //         */
 //        double Delassus[dim_M][dim_M];
 //        for (int i=0; i<dim_M; ++i){
 //            for (int j=0; j<dim_M; ++j){
 //                Delassus[i][j] = BMB(i,j);
 //            }
 //        }

 //        if (G_exist) {
 //        	const H5std_string name_data_pre = std::to_string(Group_size);
 //        	DataSet* dataset_pre = new DataSet(Matrix_G->openDataSet( name_data_pre , PredType::NATIVE_DOUBLE ));

	// 		int dim_out[2] = dataset_pre->getInMemDataSize();
 //        	double data_out[dim_out[0]][dim_out[1]];

 //        	dataset_pre->read( data_out, PredType::NATIVE_DOUBLE );

 //        	/*
 //        	 * Check if matrix already exists? (A large number of attempt to solve LCP)
 //        	 */
 //        	MatrixXd Diff( dim_M , dim_M );
 //        	for (int i=0; i<dim_M; ++i){
 //        		for (int j=0; j<dim_M; ++j){
 //        			const double val_rel = std::min( std::abs(Delassus[i][j]) , std::abs(data_out[i][j]) );
 //        			const double val_rel_a = (Delassus[i][j] - data_out[i][j])/val_rel;

 //        			Diff(i,j) = std::max( std::abs( val_rel_a ) , 0);
 //        		}
 //        	}
 //        	bool is_same_LCP = Diff.norm() < 1e-7;
 //        }

 //        /*
 //         * Create dataspace for the dataset in the file.
 //         */
 //        if (!is_same_LCP) {
	//         hsize_t dim_space_M[2];
	//         dim_space_M[0] = dim_M;
	//         dim_space_M[1] = dim_M;
	//         DataSpace fspace_M( 2, dim_space_M );
	//         /*
	//          * Create dataset and write it into the file.
	//          */
	//         DataSet* dataset_M = new DataSet(Matrix_G->createDataSet(name_matrix
	//             , PredType::NATIVE_DOUBLE, fspace_M));

	//         dataset_M->write(Delassus, PredType::NATIVE_DOUBLE);

	//         /*
	//          * Close the dataset and the file.
	//          */        
	//         delete dataset_M;

	//         /*-----------------------------------------------------------------------------------------
	//          * new dataset for relative velocities
	//          *---------------------------------------------------------------------------------------*/
	//         const H5std_string name_vector = std::to_string(Group_size+1);
	//         /*
	//          * Create dataspace for the dataset in the file.
	//          */

	//         hsize_t dim_space_V[1];
	//         dim_space_V[0] = dim_M;

	//         DataSpace fspace_V( 1, dim_space_V );
	//         /*
	//          * Create dataset and write it into the file.
	//          */

	//         DataSet* dataset_V = new DataSet(Vector_G->createDataSet(name_vector
	//             , PredType::NATIVE_DOUBLE, fspace_V));

	//         /*
	//          * Conversion Eigen -> DOUBLE
	//          */
	//         double rel_vel[dim_M];
	//         for (int i=0; i<dim_M; ++i){
	//             rel_vel[i] = BMB(i,1);
	//         }

	//         dataset_V->write(rel_vel, PredType::NATIVE_DOUBLE);

	//         /*
	//          * Close the dataset and the file.
	//          */        
	//         delete dataset_V;
 //    	}

 //        delete file;
 //   }  // end of try block
 //   // catch failure caused by the H5File operations
 //   catch( FileIException error )
 //   {
 //    error.printError();
 //   }
 //   // catch failure caused by the DataSet operations
 //   catch( DataSetIException error )
 //   {
 //    error.printError();
 //   }
 //   // catch failure caused by the DataSpace operations
 //   catch( DataSpaceIException error )
 //   {
 //    error.printError();
 //   }

   return 0;

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
// 	    error.printError();
// 	    return -1;
//    	}
//    	// catch failure caused by the DataSet operations
//    	catch( DataSetIException error )
//    	{
//     	error.printError();
//     	return -1;
//    	}
//    	// catch failure caused by the DataSpace operations
//    	catch( DataSpaceIException error )
//    	{
//     	error.printError();
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