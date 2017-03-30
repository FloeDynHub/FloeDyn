/*!
 * \file io/hdf5_floe_group_import.hpp
 * \brief Import floe group from hdf5 file
 * \author Quentin Jouet
 */

#ifndef FLOE_IO_HDF5_FLOE_GROUP_IMPORT_HPP
#define FLOE_IO_HDF5_FLOE_GROUP_IMPORT_HPP

#include <iostream>
#include "floe/generator/mesh_generator.hpp"
#include "boost/multi_array.hpp"

#include "H5Cpp.h"
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif


namespace floe { namespace io
{


template <typename TFloeGroup>
void import_floes_from_hdf5(H5std_string filename, TFloeGroup& floe_group)
{
    // using floe_type = typename TFloeList::real_type;
    using floe_type = typename TFloeGroup::floe_type;
    using geometry_type = typename floe_type::geometry_type;
    using mesh_type = typename floe_type::mesh_type;
    using point_type = typename floe_type::point_type;
    using real_type = typename floe_type::real_type;
    using static_floe_type = typename floe_type::static_floe_type;
    using array_2d_type = boost::multi_array<real_type, 2>;

    auto& floe_list = floe_group.get_floes();
    /*
     * Open the specified file and the specified dataset in the file.
     */
    H5File file( filename, H5F_ACC_RDONLY );


    // Import states
    DataSet states_dataset = file.openDataSet( "floe_states" );
    /*
    * Get dataspace of the dataset.
    */
    DataSpace states_dataspace = states_dataset.getSpace();
    /*
    * Get the number of dimensions in the dataspace.
    */
    const int rank = 3;
    /*
    * Get the dimension size of each dimension in the dataspace and
    * display them.
    */
    hsize_t states_dims_out[rank];
    states_dataspace.getSimpleExtentDims( states_dims_out, NULL);
    /*
    * Define hyperslab in the dataset; implicitly giving strike and
    * block NULL.
    */
    hsize_t i = 0;
    hsize_t      states_offset[rank] = {i, 0, 0};  // hyperslab offset in the file
    hsize_t      states_count[rank] = {1, states_dims_out[1], states_dims_out[2]};    // size of the hyperslab in the file
    states_dataspace.selectHyperslab( H5S_SELECT_SET, states_count, states_offset );
    /*
    * Define the memory dataspace.
    */
    hsize_t     states_dimsm[2] {states_dims_out[1], states_dims_out[2]};              /* memory space dimensions */
    DataSpace states_memspace( 2, states_dimsm );
    /*
    * Read data from hyperslab in the file into the hyperslab in
    * memory and display the data.
    */
    
    array_2d_type states_data_out(boost::extents[states_dims_out[1]][states_dims_out[2]]);
    states_dataset.read( states_data_out.data(), PredType::NATIVE_DOUBLE, states_memspace, states_dataspace );


    /* read shapes */
    Group floe_shape_group = file.openGroup("floe_shapes");

    // Resize floe_list
    std::size_t nb_floes = states_dims_out[1];
    std::size_t previous_nb_floes = floe_list.size();
    floe_list.resize(previous_nb_floes + nb_floes);

    for (std::size_t floe_id = 0; floe_id < nb_floes; ++floe_id)
    {
        try
        {  
            Exception::dontPrint();
            // to determine if the dataset exists in the group
            auto dataset = floe_shape_group.openDataSet(std::to_string(floe_id));

            DataSpace dataspace = dataset.getSpace();
            const int rank = 2;
            hsize_t dims_out[rank];
            dataspace.getSimpleExtentDims( dims_out, NULL);
            DataSpace memspace( 2, dims_out );
            // real_type data_out[dims_out[0]][dims_out[1]];
            array_2d_type data_out(boost::extents[dims_out[0]][dims_out[1]]);
            dataset.read( data_out.data(), PredType::NATIVE_DOUBLE, memspace, dataspace );

            // create geometry
            geometry_type shape;
            auto& boundary = shape.outer();
            for (std::size_t j = 0; j < dims_out[0]; j++)
            {
                boundary.push_back(point_type{data_out[j][0], data_out[j][1]});
            }
            // create mesh
            auto mesh = floe::generator::generate_mesh_for_shape<geometry_type, mesh_type>(shape);
            floe_type& floe = floe_list[previous_nb_floes + floe_id];
            // link static floe
            floe.attach_static_floe_ptr(std::unique_ptr<static_floe_type>(new static_floe_type()));
            auto& static_floe = floe.static_floe();
            // Attach boundary
            std::unique_ptr<geometry_type> geometry(new geometry_type(shape));
            static_floe.attach_geometry_ptr(std::move(geometry));
            // Attach mesh
            mesh_type& floe_mesh = floe.get_floe_h().m_static_mesh;
            floe_mesh = mesh;
            floe.static_floe().attach_mesh_ptr(&floe_mesh);

            floe_group.get_floe_group_h().add_floe(floe.get_floe_h());

            floe.set_state({
                {states_data_out[floe_id][0], states_data_out[floe_id][1]}, states_data_out[floe_id][2],
                {states_data_out[floe_id][3], states_data_out[floe_id][4]}, states_data_out[floe_id][5],
                {0,0}
            });
        }
        catch( GroupIException not_found_error ) {
            std::cout << "Erreur d'importation" << std::endl;
            break;
        }
    }

    // Import states
    DataSet window_dataset = file.openDataSet( "window" );
    real_type win_data[4];
    window_dataset.read( win_data, PredType::NATIVE_DOUBLE );
    floe_group.set_initial_window({{win_data[0], win_data[1], win_data[2], win_data[3]}});

    std::cout << nb_floes << " Floes, concentration : " << (int)(floe_group.initial_concentration() * 100) << "%" << std::endl;
};


}} // namespace floe::io


#endif // FLOE_IO_HDF5_FLOE_GROUP_IMPORT_HPP
