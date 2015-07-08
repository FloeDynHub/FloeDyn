/*!
 * \file variable/hdf5_writer.hpp
 * \brief HDF5 writer for output
 * \author Quentin Jouet
 */

#ifndef FLOE_IO_HDF5_WRITER_HPP
#define FLOE_IO_HDF5_WRITER_HPP

#include <iostream>
#include <vector>
#include "floe/variable/floe_group.hpp"

#include "H5Cpp.h"
#include <algorithm>
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif


namespace floe { namespace io
{

/*! HDF5Writer
 *
 * Handles floe states, shapes and time output
 *
 */


template <
    typename TFloeGroup
>
class HDF5Writer
{

public:
    template<typename T>
    using vector = std::vector<T>;
    using floe_group_type = TFloeGroup;
    using value_type = typename TFloeGroup::value_type;

    //! Default constructor.
    HDF5Writer() : m_out_file{nullptr}, m_step_count{0}, m_chunk_step_count{0}, m_flush_max_step{100}
    {
        m_data_chunk_time.reserve(m_flush_max_step);
    }

    //! Destructor
    ~HDF5Writer()
    {
        write_chunk();
        delete m_out_file;
    }

    void save_step(value_type time, const floe_group_type& floe_group);

    void write_chunk();

private:

    H5File* m_out_file;
    hsize_t m_step_count;
    hsize_t m_chunk_step_count;
    const hsize_t m_flush_max_step;
    vector<vector<vector<vector<value_type>>>> m_data_chunk_boundaries;
    vector<value_type> m_data_chunk_time;

    void write_boundaries();
    void write_time();

};


template <
    typename TFloe
>
void HDF5Writer<TFloe>::save_step(value_type time, const floe_group_type& floe_group)
{
    if (m_data_chunk_boundaries.size() == 0)
    {   
        m_data_chunk_boundaries.resize(floe_group.get_floes().size());
        for (std::size_t i = 0; i != m_data_chunk_boundaries.size(); ++i)
        {
            m_data_chunk_boundaries[i].reserve(m_flush_max_step);
        }
    }

    std::size_t floe_id = 0;
    for (auto const& floe : floe_group.get_floes())
    {
        vector<vector<value_type>> floe_step_data;
        for (auto const& pt : floe.geometry().outer())
            floe_step_data.push_back({pt.x, pt.y});
        m_data_chunk_boundaries[floe_id].push_back(floe_step_data);
        floe_id++;
    }  
    
    m_data_chunk_time[m_chunk_step_count] = time;

    m_step_count++;
    m_chunk_step_count++;

    if (m_step_count % m_flush_max_step == 0)
    {
        write_chunk();
        m_chunk_step_count = 0;
    }
};


template <typename TFloe>
void HDF5Writer<TFloe>::write_chunk() {
    try
    {   
        /*
         * Turn off the auto-printing when failure occurs so that we can
         * handle the errors appropriately
         */
        Exception::dontPrint();
        /*
         * Create a new file using H5F_ACC_TRUNC access,
         * default file creation properties, and default file
         * access properties.
         */
        if (m_out_file == nullptr)
        {
            const H5std_string  FILE_NAME( "out/out.h5" );
            m_out_file = new H5File( FILE_NAME, H5F_ACC_TRUNC );
        }

        write_boundaries();
        write_time();

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
void HDF5Writer<TFloe>::write_boundaries() {
    
    H5File& file( *m_out_file );
    const int   SPACE_DIM = 2;

    Group floe_state_group;
    try {
        floe_state_group = file.openGroup("states");
    } catch (H5::Exception& e) {
        floe_state_group = file.createGroup(H5std_string{"states"});
    }
    
    for (int i = 0; i!= m_data_chunk_boundaries.size(); ++i)
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
void HDF5Writer<TFloe>::write_time() {
    
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
    value_type data_time[m_chunk_step_count];
    for (std::size_t i = 0; i != m_chunk_step_count; ++i)
        data_time[i] = m_data_chunk_time[i];
    // Write data to the extended portion of the dataset.
    time_dataset.write(data_time, PredType::NATIVE_DOUBLE, memspace, filespace);

    // clearing buffer
    m_data_chunk_time.clear();

};


// template <
//     typename TFloe
// >
// void HDF5Writer<TFloe>::write_shapes() {
    
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


#endif // FLOE_IO_HDF5_WRITER_HPP
