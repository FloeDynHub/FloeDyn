/*!
 * \file variable/hdf5_manager.hpp
 * \brief HDF5 manager for io
 * \author Quentin Jouet
 */

#ifndef FLOE_IO_HDF5_MANAGER_HPP
#define FLOE_IO_HDF5_MANAGER_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include "floe/variable/floe_group.hpp"

#include "H5Cpp.h"
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif


namespace floe { namespace io
{

/*! HDF5Manager
 *
 * Handles floe outlines, floe states and time output
 *
 */

std::string gen_random(const int len) {
    static const char alphanum[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
        
    std::string s;

    std::srand(std::time(0));

    for (int i = 0; i < len; ++i) {
        s += alphanum[rand() % (sizeof(alphanum) - 1)];
    }

    return s;
}

template <
    typename TFloeGroup
>
class HDF5Manager
{

public:
    template<typename T>
    using vector = std::vector<T>;
    using floe_group_type = TFloeGroup;
    using value_type = typename TFloeGroup::value_type;
    using saved_state_type = std::array<value_type, 6>;

    //! Default constructor.
    HDF5Manager() :
        m_out_file_name{"io/out" + gen_random(5) + ".h5"},
        m_out_file{nullptr}, m_step_count{0}, m_chunk_step_count{0}, m_flush_max_step{100},
        m_data_chunk_time{new value_type[m_flush_max_step]}
    {
        m_data_chunk_mass_center.reserve(m_flush_max_step);
    }

    //! Destructor
    ~HDF5Manager()
    {
        if (m_chunk_step_count != 0)
            write_chunk();
        std::cout << "OUT FILE : " << m_out_file_name << std::endl;
    }

    void save_step(value_type time, const floe_group_type& floe_group);

    double recover_states(H5std_string filename, value_type time, floe_group_type& floe_group);

private:

    const std::string m_out_file_name;
    H5File* m_out_file;
    hsize_t m_step_count;
    hsize_t m_chunk_step_count;
    const hsize_t m_flush_max_step;
    vector<vector<vector<vector<value_type>>>> m_data_chunk_boundaries;
    vector<vector<saved_state_type>> m_data_chunk_states;
    // vector<value_type> m_data_chunk_time;
    value_type* m_data_chunk_time;
    vector<std::array<value_type, 2>> m_data_chunk_mass_center;

    void write_chunk();
    void write_boundaries();
    void write_states();
    void write_time();
    void write_mass_center();

};


template <
    typename TFloe
>
void HDF5Manager<TFloe>::save_step(value_type time, const floe_group_type& floe_group)
{
    auto const& floe_list = floe_group.get_floes();

    if (m_data_chunk_boundaries.size() == 0)
    {   
        m_data_chunk_boundaries.resize(floe_list.size());
        // for (std::size_t i = 0; i != m_data_chunk_boundaries.size(); ++i)
        // {
        //     m_data_chunk_boundaries[i].reserve(m_flush_max_step);
        // }
    }

    // save boundaries
    std::size_t floe_id = 0;
    for (auto const& floe : floe_list)
    {
        vector<vector<value_type>> floe_step_data;
        for (auto const& pt : floe.geometry().outer())
            floe_step_data.push_back({pt.x + floe.state().trans.x, pt.y + floe.state().trans.y});
        m_data_chunk_boundaries[floe_id].push_back(floe_step_data);
        floe_id++;
    }

    // save states
    vector<saved_state_type> states_step_data(floe_list.size());
    floe_id = 0;
    for(auto const& floe : floe_list)
    {
        states_step_data[floe_id] = {{
            floe.state().real_position().x,
            floe.state().real_position().y,
            floe.state().theta,
            floe.state().speed.x,
            floe.state().speed.y,
            floe.state().rot
        }};
        floe_id++;
    }
    m_data_chunk_states.push_back(states_step_data);

    // save time
    m_data_chunk_time[m_chunk_step_count] = time;

    // save mass center
    auto mass_center = floe_group.mass_center();
    m_data_chunk_mass_center[m_chunk_step_count] = {{mass_center.x, mass_center.y}};


    m_step_count++;
    m_chunk_step_count++;

    if (m_step_count % m_flush_max_step == 0)
    {
        write_chunk();
        m_chunk_step_count = 0;
    }
};


template <typename TFloe>
void HDF5Manager<TFloe>::write_chunk() {
    try
    {   
        /*
         * Turn off the auto-printing when failure occurs so that we can
         * handle the errors appropriately
         */
        Exception::dontPrint();

        // if (m_out_file == nullptr){}
        const H5std_string  FILE_NAME( m_out_file_name );

        if (m_step_count <= m_flush_max_step)
        {
            /*
             * Create a new file using H5F_ACC_TRUNC access,
             * default file creation properties, and default file
             * access properties.
             */
            m_out_file = new H5File( FILE_NAME, H5F_ACC_TRUNC );
            // m_out_file = H5Fcreate( FILE_NAME, H5F_ACC_TRUNC);
        } else {
            /*
             * Open the file with read/write access.
             */
            m_out_file = new H5File( FILE_NAME, H5F_ACC_RDWR );
        }

        write_boundaries();
        write_states();
        write_time();
        write_mass_center();

        delete m_out_file;

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
    // catch failure caused by the DataType operations
    catch( DataTypeIException error )
    {
        error.printError();
    }
};


template <
    typename TFloe
>
void HDF5Manager<TFloe>::write_boundaries() {
    
    H5File& file( *m_out_file );
    const int   SPACE_DIM = 2;

    Group floe_state_group;
    try {
        floe_state_group = file.openGroup("floe_outlines");
    } catch (H5::Exception& e) {
        floe_state_group = file.createGroup(H5std_string{"floe_outlines"});
    }
    
    for (std::size_t i = 0; i!= m_data_chunk_boundaries.size(); ++i)
    {
        const int   RANK = 3;
        auto& floe_chunk = m_data_chunk_boundaries[i];
        hsize_t     dimsf[RANK] = {m_step_count - m_chunk_step_count, floe_chunk[0].size(), SPACE_DIM};
        hsize_t chunk_dims[RANK] = {m_chunk_step_count, dimsf[1], dimsf[2]};

        DataSet dataset;
        try {
            dataset = floe_state_group.openDataSet(H5std_string{std::to_string(i)});
        } catch (H5::Exception& e) {
            FloatType datatype( PredType::NATIVE_DOUBLE );
            datatype.setOrder( H5T_ORDER_LE );
            hsize_t maxdims[RANK] = {H5S_UNLIMITED, dimsf[1], dimsf[2]};
            DataSpace dataspace( RANK, dimsf, maxdims );
            // Modify dataset creation property to enable chunking
            DSetCreatPropList prop;
            prop.setChunk(RANK, chunk_dims);

            dataset = floe_state_group.createDataSet(H5std_string{std::to_string(i)},datatype, dataspace, prop);
        }

        // Extend the dataset.
        dimsf[0] += m_chunk_step_count;
        dataset.extend(dimsf); 

        DataSpace filespace = dataset.getSpace();
        hsize_t offset[RANK] = {m_step_count - m_chunk_step_count, 0, 0}; // for dataset extension
        filespace.selectHyperslab(H5S_SELECT_SET, chunk_dims, offset);

        // Define memory space.
        DataSpace memspace{RANK, chunk_dims, NULL};

        value_type data[m_chunk_step_count][dimsf[1]][dimsf[2]];
        for (std::size_t j = 0; j != floe_chunk.size(); ++j)
        {
            for (std::size_t k = 0; k != floe_chunk[j].size(); ++k)
            {
                data[j][k][0] = floe_chunk[j][k][0];
                data[j][k][1] = floe_chunk[j][k][1];
            }
        }

        // Write data to the extended portion of the dataset.
        dataset.write(data, PredType::NATIVE_DOUBLE, memspace, filespace);
    }

    // clearing buffer
    m_data_chunk_boundaries.clear();

};

template <
    typename TFloe
>
void HDF5Manager<TFloe>::write_states() {
    
    H5File& file( *m_out_file );
    const int   RANK = 3;

    /* saving time */
    DataSet states_dataset;
    const hsize_t nb_floes = m_data_chunk_states[0].size();
    hsize_t     dims[RANK] = {m_step_count - m_chunk_step_count, nb_floes, 6};
    const hsize_t     chunk_dims[RANK] = {m_chunk_step_count, dims[1], dims[2]};
    try {
        states_dataset = file.openDataSet("floe_states");
    } catch (H5::Exception& e) {
        FloatType datatype( PredType::NATIVE_DOUBLE );
        datatype.setOrder( H5T_ORDER_LE );
        hsize_t maxdims[RANK] = {H5S_UNLIMITED, dims[1], dims[2]}; 
        DataSpace dataspace( RANK, dims, maxdims );
        // Modify dataset creation property to enable chunking
        DSetCreatPropList prop;
        prop.setChunk(RANK, chunk_dims);

        states_dataset = file.createDataSet(H5std_string{"floe_states"}, datatype, dataspace, prop);
    }
    // Extend the dataset.
    dims[0] += chunk_dims[0];
    states_dataset.extend(dims); 

    DataSpace filespace = states_dataset.getSpace();
    hsize_t offset[RANK] = {m_step_count - m_chunk_step_count, 0, 0};
    filespace.selectHyperslab(H5S_SELECT_SET, chunk_dims, offset);
    // Define memory space.
    DataSpace memspace{RANK, chunk_dims, NULL};
    value_type data[chunk_dims[0]][chunk_dims[1]][chunk_dims[2]];
    for (std::size_t j = 0; j != chunk_dims[0]; ++j)
    {
        for (std::size_t k = 0; k != chunk_dims[1]; ++k)
        {
            for (std::size_t l = 0; l!=6; ++l)
                data[j][k][l] = m_data_chunk_states[j][k][l];
        }
    }
    // Write data to the extended portion of the dataset.
    states_dataset.write(data, PredType::NATIVE_DOUBLE, memspace, filespace);

    // clearing buffer
    m_data_chunk_states.clear();

};


template <
    typename TFloe
>
void HDF5Manager<TFloe>::write_time() {
    
    H5File& file( *m_out_file );

    /* saving time */
    DataSet time_dataset;
    hsize_t     dimst[1] = {m_step_count - m_chunk_step_count};
    hsize_t     chunk_dimst[1] = {m_chunk_step_count};
    try {
        time_dataset = file.openDataSet("time");
    } catch (H5::Exception& e) {
        FloatType datatype( PredType::NATIVE_DOUBLE );
        datatype.setOrder( H5T_ORDER_LE );
        hsize_t maxdims[1] = {H5S_UNLIMITED}; 
        DataSpace dataspace( 1, dimst, maxdims );
        // Modify dataset creation property to enable chunking
        DSetCreatPropList prop;
        prop.setChunk(1, chunk_dimst);

        time_dataset = file.createDataSet(H5std_string{"time"}, datatype, dataspace, prop);
    }
    // Extend the dataset.
    dimst[0] += chunk_dimst[0];
    time_dataset.extend(dimst); 

    DataSpace filespace = time_dataset.getSpace();
    hsize_t offset[1] = {m_step_count - m_chunk_step_count};
    filespace.selectHyperslab(H5S_SELECT_SET, chunk_dimst, offset);
    // Define memory space.
    DataSpace memspace{1, chunk_dimst, NULL};
    time_dataset.write(m_data_chunk_time, PredType::NATIVE_DOUBLE, memspace, filespace);

};

template <
    typename TFloe
>
void HDF5Manager<TFloe>::write_mass_center() {
    
    H5File& file( *m_out_file );
    const int   RANK = 2;

    /* saving mass center */
    DataSet dataset;
    hsize_t     dimst[RANK] = {m_step_count - m_chunk_step_count, 2};
    hsize_t     chunk_dims[RANK] = {m_chunk_step_count, 2};
    try {
        dataset = file.openDataSet("mass_center");
    } catch (H5::Exception& e) {
        FloatType datatype( PredType::NATIVE_DOUBLE );
        datatype.setOrder( H5T_ORDER_LE );
        hsize_t maxdims[RANK] = {H5S_UNLIMITED, 2}; 
        DataSpace dataspace( RANK, dimst, maxdims );
        // Modify dataset creation property to enable chunking
        DSetCreatPropList prop;
        prop.setChunk(RANK, chunk_dims);

        dataset = file.createDataSet(H5std_string{"mass_center"}, datatype, dataspace, prop);
    }
    // Extend the dataset.
    dimst[0] += chunk_dims[0];
    dataset.extend(dimst); 

    DataSpace filespace = dataset.getSpace();
    hsize_t offset[RANK] = {m_step_count - m_chunk_step_count, 0};
    filespace.selectHyperslab(H5S_SELECT_SET, chunk_dims, offset);
    // Define memory space.
    DataSpace memspace{RANK, chunk_dims, NULL};

    value_type data[chunk_dims[0]][chunk_dims[1]];
    for (std::size_t j = 0; j != chunk_dims[0]; ++j)
    {
        for (std::size_t k = 0; k != chunk_dims[1]; ++k)
        {
            data[j][k] = m_data_chunk_mass_center[j][k];
        }
    }

    dataset.write(data, PredType::NATIVE_DOUBLE, memspace, filespace);

};

template <
    typename TFloe
>
double HDF5Manager<TFloe>::recover_states(H5std_string filename, value_type time, floe_group_type& floe_group) {
    
    /*
     * Open the specified file and the specified dataset in the file.
     */
    H5File file( filename, H5F_ACC_RDONLY );

    DataSet time_dataset = file.openDataSet( "time" );
    /*
    * Get dataspace of the dataset.
    */
    DataSpace time_dataspace = time_dataset.getSpace();
    /*
    * Get the dimension size of each dimension in the dataspace and
    * display them.
    */
    hsize_t dims_out[1];
    time_dataspace.getSimpleExtentDims( dims_out, NULL);
    value_type data_time[dims_out[0]];
    for (std::size_t j = 0; j!= dims_out[0]; ++j)
        data_time[j] = 0;
     /*
    * Define the memory dataspace.
    */
    DataSpace time_memspace( 1, dims_out );
    /*
    * Read data from the file
    */
    time_dataset.read( data_time, PredType::NATIVE_DOUBLE, time_memspace, time_dataspace );

    hsize_t i = 0;
    while (data_time[i] < time && i < dims_out[0])
        ++i;
    --i;


    {
    DataSet dataset = file.openDataSet( "floe_states" );
    /*
    * Get dataspace of the dataset.
    */
    DataSpace dataspace = dataset.getSpace();
    /*
    * Get the number of dimensions in the dataspace.
    */
    const int rank = 3;
    /*
    * Get the dimension size of each dimension in the dataspace and
    * display them.
    */
    hsize_t dims_out[rank];
    dataspace.getSimpleExtentDims( dims_out, NULL);
    /*
    * Define hyperslab in the dataset; implicitly giving strike and
    * block NULL.
    */
    hsize_t      offset[rank] = {i, 0, 0};  // hyperslab offset in the file
    hsize_t      count[rank] = {1, dims_out[1], dims_out[2]};    // size of the hyperslab in the file
    dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );
    /*
    * Define the memory dataspace.
    */
    hsize_t     dimsm[2] {dims_out[1], dims_out[2]};              /* memory space dimensions */
    DataSpace memspace( 2, dimsm );
    /*
    * Read data from hyperslab in the file into the hyperslab in
    * memory and display the data.
    */
    value_type data_out[dims_out[1]][dims_out[2]];
    dataset.read( data_out, PredType::NATIVE_DOUBLE, memspace, dataspace );

    std::size_t floe_id = 0;
    for (auto& floe : floe_group.get_floes())
    {
        floe.set_state({
            {data_out[floe_id][0], data_out[floe_id][1]}, data_out[floe_id][2],
            {data_out[floe_id][3], data_out[floe_id][4]}, data_out[floe_id][5]
        });
        floe_id++;
    }

    return data_time[i];

    }

};


// template <
//     typename TFloe
// >
// void HDF5Manager<TFloe>::write_shapes() {
    
//     H5File& file( *m_out_file );

//     /* write shapes */
//     try {
//         Group floe_shape_group = file.openGroup("shapes");
//     } catch (H5::Exception& e) {
//         /* Create group for floe shapes */
//         Group floe_shape_group(file.createGroup(H5std_string{"shapes"}));

//         const int   RANK = 2;
//         FloatType datatype( PredType::NATIVE_DOUBLE );
//         datatype.setOrder( H5T_ORDER_LE );
//         hsize_t     dimsf[2];              // dataset dimensions
//         dimsf[1] = SPACE_DIM;
//         for (std::size_t i=0; i!=floe_group.get_floes().size(); ++i)
//         {
//             auto& boundary = floe_group.get_floes()[i].geometry().outer();
//             dimsf[0] = boundary.size();
//             DataSpace dataspace( RANK, dimsf );
//             /*
//              * Create a nFew dataset within the file using defined dataspace and
//              * datatype and default dataset creation properties.
//              */
//             DataSet dataset = floe_shape_group.createDataSet(H5std_string{std::to_string(i)},datatype, dataspace);
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
//     }

// };


}} // namespace floe::io


#endif // FLOE_IO_HDF5_MANAGER_HPP
