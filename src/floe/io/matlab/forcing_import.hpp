/*!
 * \file floe/io/matlab/forcing_import.hpp
 * \brief Import of ocean/wind speeds from an external Matlab file.
 * \author Samuel Brenner
 * \date 07/30/2022 (updated 12/07/2022 to include ug,vg)
 * \ adapted from floe/io/matlab/topaz_import.hpp by Quentin Jouet
 */

#ifndef FLOE_IO_MATLAB_FORCING_IMPORT_HPP
#define FLOE_IO_MATLAB_FORCING_IMPORT_HPP

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
template <typename T, typename Index>
void read_forcing_from_file( std::string const& file_name,
                             const T* &x, const T* &y, const T* &t,
                             const T* &u, const T* &v,
                             const T* &ug, const T* &vg,                             
                             Index &nd)
{
    // std::string file_name = "io/inputs/input_forcing.mat";
    std::cout << "Reading " << file_name.c_str() <<" ...\n" ;
    mat_t *matfp;
    matfp = Mat_Open( file_name.c_str(), MAT_ACC_RDONLY );
    if ( matfp == nullptr )
    {
        throw std::ios_base::failure("Error opening MAT file");
    }
    // READ MATLAB VARIABLES
    // note: if variable missing from .mat file, it will be read in here as NULL
    matvar_t *xmat = Mat_VarRead(matfp,"x");
    matvar_t *ymat = Mat_VarRead(matfp,"y");
    matvar_t *tmat = Mat_VarRead(matfp,"t");
    matvar_t *umat = Mat_VarRead(matfp,"u");
    matvar_t *vmat = Mat_VarRead(matfp,"v");    
    matvar_t *ugmat = Mat_VarRead(matfp,"ug");
    matvar_t *vgmat = Mat_VarRead(matfp,"vg");  
    // matvar_t *uamat = Mat_VarRead(matfp,"ua");
    // matvar_t *vamat = Mat_VarRead(matfp,"va");  


    // ASSIGN DATA READ
    x = static_cast<const double*>(xmat->data);
    y = static_cast<const double*>(ymat->data);
    t = static_cast<const double*>(tmat->data);
    if (umat!=NULL && vmat!=NULL){
        u = static_cast<const double*>(umat->data);
        v = static_cast<const double*>(vmat->data);
        ug = u;
        vg = v;  
    }
    if (ugmat!=NULL && vgmat!=NULL){
        ug = static_cast<const double*>(ugmat->data);
        vg = static_cast<const double*>(vgmat->data);
    }
    // if (uamat!=NULL && vamat!=NULL){
    //     ua = static_cast<const double*>(uamat->data);
    //     va = static_cast<const double*>(vamat->data);
    // }

    Mat_Close(matfp);

    // Get input file dimensions
    int N = (xmat->dims[1]);
    int M = (ymat->dims[1]);
    int L = (tmat->dims[1]);
    nd[0] = N; nd[1] = M; nd[2] = L; 


}


}}} // namespace floe::io::matlab

#endif // FLOE_IO_MATLAB_FORCING_IMPORT_HPP

