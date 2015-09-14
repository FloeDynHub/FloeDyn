/*!
 * \file floe/io/matlab/list_so_import.hpp
 * \brief Import of list_so structure for a Matlab file.
 * \author Roland Denis
 */

#ifndef FLOE_IO_MATLAB_LIST_SO_IMPORT_HPP
#define FLOE_IO_MATLAB_LIST_SO_IMPORT_HPP

#include <ios>
#include <string>
#include <vector>
#include <matio.h>
#include "floe/io/matlab/list_so.hpp"

namespace floe { namespace io { namespace matlab
{

using namespace std;

template <typename T> void read_list_so( matvar_t* var, MatlabListSolid<T>& list_so );
template <typename T> void read_geo( matvar_t* var, std::vector< MatlabSolid<T> >& geo );
template <typename T> void read_mov( matvar_t* var, std::vector< MatlabMovement<T> >& mov );

template <typename T> void read_value( matvar_t* var, T& value );
template <typename T> void read_value_list( matvar_t* var, std::vector<T>& value_list );
template <typename T> void read_point( matvar_t* var, point_type<T>& value);
template <typename T> void read_point_list( matvar_t* var, std::vector< point_type<T> >& point_list );
template <typename T> void read_border( matvar_t* var, border_type<T>& border );
template <typename T> void read_mark( matvar_t* var, mark_type<T>& mark );
template <typename T> void read_mesh( matvar_t* var, mesh_type<T>& mesh );
template <typename T> void read_triangle_list( matvar_t* var, std::vector< std::array<T, 3> >& triangle_list );
template <typename T> void read_own_disk( matvar_t* var, own_disk_type<T>& own_disk );
template <typename T> void read_pt_own_mark( matvar_t* var, pt_own_mark_type<T>& pt_own_mark );
template <typename T> void read_disk_rep_list( matvar_t* var, std::vector< disk_rep_type<T> >& disk_rep_list );
template <typename T> void read_disks( matvar_t* var, disks_type<T>& disks );
template <typename T> void read_quadrature_list( matvar_t* var, std::vector< quadrature_type<T> >& quadrature_list );


/*! Reads list_so from a .mat file
 *
 * \tparam T Value type used by matlab.
 * \param file_name File name.
 * \param list_so   Object to store the result.
 * \param variable_name Name of the list_so variable in the Matlab file.
 */
template <typename T>
void read_list_so_from_file( std::string const& file_name, MatlabListSolid<T>& list_so, std::string const& variable_name = "list_so_str" )
{
    mat_t *matfp;
    matvar_t *matvar;

    // Opening file
    matfp = Mat_Open( file_name.c_str(), MAT_ACC_RDONLY );
    if ( matfp == nullptr )
    {
        throw std::ios_base::failure("Error opening MAT file \"" + file_name + "\"");
    }

    // Reading list_so variable 
    matvar = Mat_VarRead( matfp, variable_name.c_str() );
    read_list_so(  matvar, list_so );
    Mat_VarFree(matvar);

    Mat_Close(matfp);
}

template <typename T>
void read_list_so( matvar_t* var, MatlabListSolid<T>& list_so )
{
    if ( var == nullptr )
        throw std::ios_base::failure("list_so not found");

    read_geo( Mat_VarGetStructFieldByName(var, "geo", 0),  list_so.geo );
    read_mov( Mat_VarGetStructFieldByName(var, "mov", 0),  list_so.mov );
    read_value( Mat_VarGetStructFieldByName(var, "cinetic", 0), list_so.cinetic);
}

template <typename T>
void read_geo( matvar_t* var, std::vector< MatlabSolid<T> >& geo )
{
    if ( var == nullptr )
        throw std::ios_base::failure("list_so.geo not found");
    
    std::size_t size = var->dims[0] * var->dims[1];
    geo.resize(size);
    for ( std::size_t i = 0; i < size; ++i )
    {
        read_value(     Mat_VarGetStructFieldByName(var, "tau", i), geo[i].tau );
        read_value(     Mat_VarGetStructFieldByName(var, "size_maille", i), geo[i].size_maille );
        read_point(     Mat_VarGetStructFieldByName(var, "center_mass", i), geo[i].center_mass );
        read_mark(      Mat_VarGetStructFieldByName(var, "mark", i), geo[i].mark );
        read_border(    Mat_VarGetStructFieldByName(var, "border", i), geo[i].border );
        read_pt_own_mark(   Mat_VarGetStructFieldByName(var, "pt_own_mark", i), geo[i].pt_own_mark );
        read_point(     Mat_VarGetStructFieldByName(var, "center_disk", i), geo[i].center_disk );
        read_value(     Mat_VarGetStructFieldByName(var, "radius", i), geo[i].radius );
        read_value(     Mat_VarGetStructFieldByName(var, "area", i), geo[i].area );
        read_disks(     Mat_VarGetStructFieldByName(var, "disks", i), geo[i].disks );
        read_value(     Mat_VarGetStructFieldByName(var, "d_c", i), geo[i].d_c );
    }

}

template <typename T> 
void read_mov( matvar_t* var, std::vector< MatlabMovement<T> >& mov )
{
    if (var == nullptr )
        throw std::ios_base::failure("list_so.mov not found");

    std::size_t size = var->dims[0] * var->dims[1];
    mov.resize(size);
    for ( std::size_t i = 0; i < size; ++i )
    {
        read_point(     Mat_VarGetStructFieldByName(var, "center", i), mov[i].center );
        read_point(     Mat_VarGetStructFieldByName(var, "vcenter", i), mov[i].vcenter );
        read_value(     Mat_VarGetStructFieldByName(var, "theta", i), mov[i].theta );
        read_value(     Mat_VarGetStructFieldByName(var, "vtheta", i), mov[i].vtheta );
        read_quadrature_list( Mat_VarGetStructFieldByName(var, "pt_reel", i), mov[i].pt_reel );
        read_value_list( Mat_VarGetStructFieldByName(var, "detJ", i), mov[i].detJ );
        read_value(     Mat_VarGetStructFieldByName(var, "mass_tot", i), mov[i].mass_tot );
        read_value_list( Mat_VarGetStructFieldByName(var, "mass", i), mov[i].mass );
        read_value_list( Mat_VarGetStructFieldByName(var, "area", i), mov[i].area );
        read_value(     Mat_VarGetStructFieldByName(var, "denom", i), mov[i].denom );
        read_value(     Mat_VarGetStructFieldByName(var, "time", i), mov[i].time );
        read_value(     Mat_VarGetStructFieldByName(var, "timecurr", i), mov[i].timecurr );
        read_value(     Mat_VarGetStructFieldByName(var, "obs", i), mov[i].obs );
    }

}

template <typename T>
void read_value( matvar_t* var, T& value )
{
    if ( var == nullptr )
        throw std::ios_base::failure("value not found");

    value = *static_cast<T*>(var->data);
}

template <typename T>
void read_value_list( matvar_t* var, std::vector<T>& value_list )
{
    if ( var == nullptr )
        throw std::ios_base::failure("value list not found");

    std::size_t size = var->dims[0] * var->dims[1];

    value_list.resize(var->dims[0]);
    for ( std::size_t i = 0; i < size; ++i )
        value_list[i] = static_cast<T*>(var->data)[i];
}

template <typename T>
void read_point( matvar_t* var, point_type<T>& value)
{
    if ( var == nullptr )
        throw std::ios_base::failure("point not found");

    std::size_t size = var->dims[0]*var->dims[1];
    if (size != 2 )
        throw std::ios_base::failure("2D point has wrong dimension");

    for ( std::size_t i = 0; i < size; ++i )
        value[i] = static_cast<T*>(var->data)[i];
}

template <typename T>
void read_point_list( matvar_t* var, std::vector< point_type<T> >& point_list )
{
    if ( var == nullptr )
        throw std::ios_base::failure("point list not found");

    std::size_t size = var->dims[0] * var->dims[1];

    point_list.resize(var->dims[0]);
    for ( std::size_t i = 0; i < size; ++i )
        point_list[i%var->dims[0]][i/var->dims[0]] = static_cast<T*>(var->data)[i];
}


template <typename T>
void read_border( matvar_t* var, border_type<T>& border )
{
    if ( var == nullptr )
        throw std::ios_base::failure("border not found");

    read_point_list( Mat_VarGetStructFieldByName(var, "S", 0), border.S);
    read_point_list( Mat_VarGetStructFieldByName(var, "ssdif", 0), border.ssdif);
}

template <typename T>
void read_mark( matvar_t* var, mark_type<T>& mark )
{
    if ( var == nullptr )
        throw std::ios_base::failure("mark not found");

    if (var->dims[0] != 3 || var->dims[1] != 2 )
        throw std::ios_base::failure("mark (3x2) has wrong dimension");

    for ( std::size_t i = 0; i < 6; ++i )
        mark[i%3][i/3] = static_cast<T*>(var->data)[i];
}



template <typename T> 
void read_mesh( matvar_t* var, mesh_type<T>& mesh )
{
    if ( var == nullptr )
        throw std::ios_base::failure("mesh not found");

    read_point_list( Mat_VarGetStructFieldByName(var, "vertices", 0), mesh.vertices);
    read_point_list( Mat_VarGetStructFieldByName(var, "border", 0), mesh.border);
    read_triangle_list( Mat_VarGetStructFieldByName(var, "triangles", 0), mesh.triangles);
}

template <typename T> 
void read_triangle_list( matvar_t* var, std::vector< std::array<T, 3> >& triangle_list )
{
    if ( var == nullptr )
        throw std::ios_base::failure("triangle list not found");
    
    std::size_t size = var->dims[0] * var->dims[1];

    triangle_list.resize(var->dims[0]);
    for ( std::size_t i = 0; i < size; ++i )
        triangle_list[i%var->dims[0]][i/var->dims[0]] = static_cast<T*>(var->data)[i];
}

template <typename T> 
void read_own_disk( matvar_t* var, own_disk_type<T>& own_disk )
{
    if ( var == nullptr )
        throw std::ios_base::failure("own_disk not found");

    read_point( Mat_VarGetStructFieldByName(var, "center_disk", 0), own_disk.center_disk );
}

template <typename T>
void read_pt_own_mark( matvar_t* var, pt_own_mark_type<T>& pt_own_mark )
{
    if ( var == nullptr )
        throw std::ios_base::failure("pt_own_mark not found");

    read_mesh( Mat_VarGetStructFieldByName(var, "S", 0), pt_own_mark.S );
    read_own_disk( Mat_VarGetStructFieldByName(var, "D", 0), pt_own_mark.D );
}

template <typename T> 
void read_disk_rep_list( matvar_t* var, std::vector< disk_rep_type<T> >& disk_rep_list )
{
    if ( var == nullptr )
        throw std::ios_base::failure("disks.rep not found");

    std::size_t size = var->dims[0] * var->dims[1];

    disk_rep_list.resize(var->dims[0]);
    for (std::size_t i = 0; i < size; ++i )
    {
        const T value = static_cast<T*>(var->data)[i]; 
        if (i < var->dims[0]) 
            disk_rep_list[i%var->dims[0]].id = value;
        else
            disk_rep_list[i%var->dims[0]].radius = value;
    }
}

template <typename T> 
void read_disks( matvar_t* var, disks_type<T>& disks )
{
    if ( var == nullptr )
        throw std::ios_base::failure("disks.rep not found");

    read_value( Mat_VarGetStructFieldByName(var, "nC", 0), disks.nC );
    read_disk_rep_list( Mat_VarGetStructFieldByName(var, "rep", 0), disks.rep );
}

template <typename T> 
void read_quadrature_list( matvar_t* var, std::vector< quadrature_type<T> >& quadrature_list )
{
    if ( var == nullptr )
        throw std::ios_base::failure("pt_reel not found");
    if ( var->rank != 3 )
        throw std::ios_base::failure("pt_reel is not a 3D array");
    if ( var->dims[0] != 0 && ( var->dims[1] != 3 || var->dims[2] != 2 ) )
        throw std::ios_base::failure("pt_reel row is not a 3x2 array");

    std::size_t size = var->dims[0] * var->dims[1] * var->dims[2];

    quadrature_list.resize(var->dims[0]);
    for ( std::size_t pos = 0; pos < size; ++pos )
    {
        const auto i = pos % var->dims[0];
        const auto j = ( pos / var->dims[0] ) % var->dims[1];
        const auto k = ( pos / var->dims[0] ) / var->dims[1];
        quadrature_list[i].points[j][k] = static_cast<T*>(var->data)[pos];
    }
}

}}} // namespace floe::io::matlab

#endif // FLOE_IO_MATLAB_LIST_SO_IMPORT_HPP

