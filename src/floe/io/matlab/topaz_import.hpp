/*!
 * \file floe/io/matlab/topaz_import.hpp
 * \brief Import of ocean/wind speeds from a topaz Matlab file.
 * \author Quentin
 */

#ifndef FLOE_IO_MATLAB_TOPAZ_IMPORT_HPP
#define FLOE_IO_MATLAB_TOPAZ_IMPORT_HPP

#include <ios>
#include <string>
#include <vector>
#include <matio.h>

#include "floe/geometry/geometries/point.hpp" 

#include <iostream> // DEBUG

namespace floe { namespace io { namespace matlab
{

using namespace std;


/*! Reads ocean and air data from a .mat file
 *
 * \tparam T Value type used by matlab.
 * \param file_name File name.
 * \param ocean_data reference to vector that will store ocean data.
 * \param air_data reference to vector that will store ocean data.
 */
template <typename TPoint>
void read_topaz_from_file( std::string const& file_name,
                           std::vector<TPoint>& ocean_data,
                           std::vector<TPoint>& air_data)
{
    using T = decltype(TPoint::x); // Coordinate type

    mat_t *matfp;

    // Opening file
    matfp = Mat_Open( file_name.c_str(), MAT_ACC_RDONLY );
    if ( matfp == nullptr )
    {
        throw std::ios_base::failure("Error opening MAT file \"" + file_name + "\"");
    }

    matvar_t *x_air;
    matvar_t *y_air;
    matvar_t *x_ocean;
    matvar_t *y_ocean;

    x_air = Mat_VarRead(matfp,"var03");
    y_air = Mat_VarRead(matfp,"var04");
    x_ocean = Mat_VarRead(matfp,"var05");
    y_ocean = Mat_VarRead(matfp,"var06");

    // checking dimensions
    for (auto* matvar : {x_air, y_air, x_ocean, y_ocean})
    {
        if ( NULL == matvar ) {
            fprintf(stderr,"Variable not found, or error reading MAT file\n");
        } else {
            if ( matvar->rank != 2 || (matvar->dims[0] > 1 && matvar->dims[1] > 1) )
                fprintf(stderr,"Variable %s is not a vector!\n", matvar->name);
        }
    }

    if ((x_air->dims[0] != y_air->dims[0]) || (x_ocean->dims[0] != y_ocean->dims[0]))
        fprintf(stderr,"x and y vector sizes are not equal !\n");

    for ( std::size_t i = 0; i < x_air->dims[0]; ++i )
        ocean_data.push_back({static_cast<T*>(x_ocean->data)[i], static_cast<T*>(y_ocean->data)[i]});

    for ( std::size_t i = 0; i < x_ocean->dims[0]; ++i )
        air_data.push_back({static_cast<T*>(x_air->data)[i], static_cast<T*>(y_air->data)[i]});

    // freeing memory
    for (auto* matvar : {x_air, y_air, x_ocean, y_ocean})
        Mat_VarFree(matvar);

    Mat_Close(matfp);
}


}}} // namespace floe::io::matlab

#endif // FLOE_IO_MATLAB_TOPAZ_IMPORT_HPP

