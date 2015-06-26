#include "../tests/catch.hpp"

#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"
#include "floe/variable/floe_group.hpp"

#include <iostream>
#include <string>
#include "H5Cpp.h"
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif
const H5std_string  FILE_NAME( "/Users/Serge/Desktop/Visu/test.h5" );
// const H5std_string  DATASET_NAME( "floe" );
// const H5std_string  DATASET2_NAME( "IntArray2" );
// const int   NX = 5;                    // dataset dimensions
const int   NY = 2;
const int   RANK = 2;

namespace ff = floe::floes;

TEST_CASE( "Test hdf5 output", "[io]" )
{

    using namespace floe::variable;
    using namespace std;
    using floe_type = ff::KinematicFloe<ff::StaticFloe<double>>;
    using point_type = typename floe_type::point_type;
    using value_type = typename floe_type::value_type;
    using mesh_type = floe_type::mesh_type;
    FloeGroup<floe_type> F;
    // Import floes from Matlab configuration
    std::string mat_file_name = "tests/floe/io/matlab/r1day_set_up_250sm_sz_60_list_so_350_str.mat";
    // std::string mat_file_name = "tests/floe/io/matlab/config_q2_str.mat";
    F.load_matlab_config(mat_file_name);

    F.out_hdf5(0);
}


// try
//     {
//         /*
//          * Turn off the auto-printing when failure occurs so that we can
//          * handle the errors appropriately
//          */
//         Exception::dontPrint();
//         /*
//          * Create a new file using H5F_ACC_TRUNC access,
//          * default file creation properties, and default file
//          * access properties.
//          */
//         H5File file( FILE_NAME, H5F_ACC_TRUNC );
//         /*
//          * Define the size of the array and create the data space for fixed
//          * size dataset.
//          */

//         FloatType datatype( PredType::NATIVE_DOUBLE );
//         datatype.setOrder( H5T_ORDER_LE );
//         hsize_t     dimsf[2];              // dataset dimensions
//         // dimsf[0] = NX;
//         dimsf[1] = NY;

//         for (std::size_t i=0; i!=F.get_floes().size(); ++i)
//         {
//             auto& boundary = F.get_floes()[i].geometry().outer();
//             dimsf[0] = boundary.size();
//             DataSpace dataspace( RANK, dimsf );
//             /*
//              * Create a new dataset within the file using defined dataspace and
//              * datatype and default dataset creation properties.
//              */
//             DataSet dataset = file.createDataSet( H5std_string{std::to_string(i)}, datatype, dataspace );
//             value_type data[dimsf[0]][dimsf[1]];
//             for (std::size_t j = 0; j!= dimsf[0]; ++j)
//             {
//                 data[j][0] = boundary[j].x;
//                 data[j][1] = boundary[j].y;
//             }
//             /*
//              * Write the data to the dataset using default memory space, file
//              * space, and transfer properties.
//              */
//             dataset.write( data, PredType::NATIVE_DOUBLE );
//         }
        
//         // Group group(file.createGroup("/MyGroup"));
//         // DataSet dataset3 = group.createDataSet(DATASET_NAME,datatype, dataspace);

//     }  // end of try block
//     // catch failure caused by the H5File operations
//     catch( FileIException error )
//     {
//         error.printError();
//         // return -1;
//     }
//     // catch failure caused by the DataSet operations
//     catch( DataSetIException error )
//     {
//         error.printError();
//         // return -1;
//     }
//     // catch failure caused by the DataSpace operations
//     catch( DataSpaceIException error )
//     {
//         error.printError();
//         // return -1;
//     }
//     // catch failure caused by the DataSpace operations
//     catch( DataTypeIException error )
//     {
//         error.printError();
//         // return -1;
//     }
