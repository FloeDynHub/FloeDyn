/*!
 * \file floe/io/matlab/pze_import.hpp
 * \brief Import box
 * \author Quentin
 */

#ifndef FLOE_IO_MATLAB_PZE_IMPORT_HPP
#define FLOE_IO_MATLAB_PZE_IMPORT_HPP

#include <ios>
#include <string>
#include <vector>
#include <array>
#include <matio.h>
#include <algorithm>


#include <iostream> // DEBUG

namespace floe { namespace io { namespace matlab
{

using namespace std;


/*! Reads ocean and air data from a .mat file
 *
 * \tparam T Value type used by matlab.
 * \param file_name File name.
 */
template<typename T=double>
std::array<T, 4> read_pze_from_file( std::string const& file_name )
{
    mat_t *matfp;

    // Opening file
    matfp = Mat_Open( file_name.c_str(), MAT_ACC_RDONLY );
    if ( matfp == nullptr )
    {
        throw std::ios_base::failure("Error opening MAT file \"" + file_name + "\"");
    }

    matvar_t *pze;

    pze = Mat_VarRead(matfp,"Pze");


    // checking dimensions
    if ( NULL == pze ) {
        fprintf(stderr,"Variable not found, or error reading MAT file\n");
        return {{0,0,0,0}};
    }

    T min_x, min_y, max_x, max_y;
    min_x = min_y = std::numeric_limits<T>::max();
    max_x = max_y = - std::numeric_limits<T>::max();

    for ( std::size_t i = 0; i < pze->dims[0] * pze->dims[1] / 2; i++ )
    {
        min_x = std::min(min_x, static_cast<T*>(pze->data)[i]);
        max_x = std::max(max_x, static_cast<T*>(pze->data)[i]);
    }
    for ( std::size_t i = pze->dims[0] * pze->dims[1] / 2; i < pze->dims[0] * pze->dims[1]; i++ )
    {
        min_y = std::min(min_y, static_cast<T*>(pze->data)[i]);
        max_y = std::max(max_y, static_cast<T*>(pze->data)[i]);
    }

    // std::cout << "***PZE***";

    // freeing memory
    Mat_VarFree(pze);

    Mat_Close(matfp);

    return {{min_x, max_x, min_y, max_y}};
}

template <typename TTopology>
TTopology topology_from_file( std::string const& file_name )
{   
    auto A = read_pze_from_file(file_name);
    return TTopology{A[0], A[1], A[2], A[3]};
}

template<typename T=double>
double ocean_window_area_from_file( std::string const& file_name )
{   
    auto A = read_pze_from_file(file_name);
    return (A[1] - A[0]) * (A[3] - A[2]);
}


}}} // namespace floe::io::matlab

#endif // FLOE_IO_MATLAB_PZE_IMPORT_HPP

