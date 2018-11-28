/*!
 * \file floes/hdf5_manager.hpp
 * \brief HDF5 manager for io
 * \author Quentin Jouet
 */

// #ifndef FLOE_IO_HDF5_MANAGER_DEF_HPP
// #define FLOE_IO_HDF5_MANAGER_DEF_HPP

#include <mpi.h>
#include "floe/io/hdf5_manager.h"
#include "floe/utils/random.hpp"

namespace floe { namespace io
{

//! Default constructor.
template <typename TFloeGroup, typename TDynamicsMgr>
HDF5Manager<TFloeGroup, TDynamicsMgr>::HDF5Manager(floe_group_type const& floe_group) :
    m_out_file_name{"io/outputs/out_" + floe::random::gen_random(5) + ".h5"},
    m_out_file{nullptr}, m_step_count{0}, m_chunk_step_count{0}, m_flush_max_step{2},
    m_floe_group{&floe_group},
    m_data_chunk_states(boost::extents[0][0][0]),
    m_data_chunk_time{new real_type[m_flush_max_step]},
    m_data_chunk_mass_center(boost::extents[m_flush_max_step][2]),
    m_data_chunk_OBL_speed(boost::extents[m_flush_max_step][2]),
    m_data_chunk_kinE{new real_type[m_flush_max_step]},
    m_out_step{0}, m_next_out_limit{0}
    {}

//! Definition of the destructor:
template <typename TFloeGroup, typename TDynamicsMgr>
HDF5Manager<TFloeGroup, TDynamicsMgr>::~HDF5Manager()
{
    flush();
    if (m_step_count) std::cout << "OUT FILE : " << m_out_file_name << std::endl;
    // if (m_data_chunk_time) {delete[] m_data_chunk_time;}
    // if (m_data_chunk_kinE) {delete[] m_data_chunk_kinE;} // BUG: sometimes freed allocated memory not allocated!!
}

template <typename TFloeGroup, typename TDynamicsMgr>
void HDF5Manager<TFloeGroup, TDynamicsMgr>::save_step_if_needed(real_type time, const dynamics_mgr_type& dynamics_manager)
{
    if (this->need_step_output(time))
    {
        this->save_step(time, dynamics_manager);
        this->update_next_out_limit();
    }
}


template <typename TFloeGroup, typename TDynamicsMgr>
void HDF5Manager<TFloeGroup, TDynamicsMgr>::save_step(real_type time, const dynamics_mgr_type& dynamics_manager)
{
    floe_group_type const& floe_group = *m_floe_group;
    // auto const& floe_list = floe_group.get_floes();
    /*
    if (m_data_chunk_boundaries.size() == 0)
    {   
        m_data_chunk_boundaries.resize(floe_list.size());
    }
    */
    // save boundaries
    // std::size_t floe_id = 0;
    /*
    for (auto const& floe : floe_list)
    {
        vector<vector<real_type>> floe_step_data;
        for (auto const& pt : floe.geometry().outer())
            floe_step_data.push_back({pt.x + floe.state().trans.x, pt.y + floe.state().trans.y});
        m_data_chunk_boundaries[floe_id].push_back(floe_step_data);
        floe_id++;
    }
    */
    // save states
    if (m_data_chunk_states.size() == 0) m_data_chunk_states.resize(boost::extents[m_flush_max_step][this->nb_considered_floes()][9]);
    for(std::size_t id = 0; id < this->nb_considered_floes(); id++)
    {
        auto const& floe = this->get_floe(id);
        int k{0};
        for (auto& val: {
            floe.state().real_position().x,
            floe.state().real_position().y,
            floe.state().theta,
            floe.state().speed.x,
            floe.state().speed.y,
            floe.state().rot,
            floe.total_received_impulse(),
            floe.state().pos.x,
            floe.state().pos.y
        }){
            m_data_chunk_states[m_chunk_step_count][id][k++] = val;
        }
    }

    // save time
    m_data_chunk_time[m_chunk_step_count] = time;

    // save mass center
    auto mass_center = floe_group.mass_center();
    m_data_chunk_mass_center[m_chunk_step_count][0] = mass_center.x;
    m_data_chunk_mass_center[m_chunk_step_count][1] = mass_center.y;

    // save OBL speed
    auto OBL_speed = dynamics_manager.OBL_speed();
    m_data_chunk_OBL_speed[m_chunk_step_count][0] = OBL_speed.x;
    m_data_chunk_OBL_speed[m_chunk_step_count][1] = OBL_speed.y;

    // save Kinetic Energy:
    m_data_chunk_kinE[m_chunk_step_count] = floe_group.kinetic_energy();

    m_step_count++;
    m_chunk_step_count++;

    if (m_step_count % m_flush_max_step == 0) // TODO change test to (m_chunk_step_count == m_flush_max_step) ?
    {
        flush();
        m_chunk_step_count = 0;
    }
};


template <typename TFloeGroup, typename TDynamicsMgr>
void HDF5Manager<TFloeGroup, TDynamicsMgr>::flush() {
    if (m_chunk_step_count == 0)
        return;
    try
    {   
        /*
         * Turn off the auto-printing when failure occurs so that we can
         * handle the errors appropriately
         */
        Exception::dontPrint();

        const H5std_string  FILE_NAME( m_out_file_name );

        if (m_step_count == m_chunk_step_count)
        {
            /*
             * Create a new file using H5F_ACC_TRUNC access,
             * default file creation properties, and default file
             * access properties.
             */
            m_out_file = new H5File( FILE_NAME.c_str(), H5F_ACC_TRUNC );
        } else {
            /*
             * Open the file with read/write access.
             */
            m_out_file = new H5File( FILE_NAME.c_str(), H5F_ACC_RDWR );
        }

        try { m_out_file->openGroup("floe_shapes"); }
        catch (...) { write_shapes(); }

        try { m_out_file->openDataSet("window"); }
        catch (...) { write_window(); }

        // write_boundaries();
        write_states();
        write_time();
        write_mass_center();
        write_OBL_speed();
        write_kinE();

        // Close the file after each flush to keep a valid ouput even if program crashes
        delete m_out_file;

    }  // end of try block
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        #ifdef MULTIOUTPUT
            char procname[MPI_MAX_PROCESSOR_NAME];
            int resultlength;
            int code = MPI_Get_processor_name(procname,&resultlength);
            int rank;
            MPI_Comm_rank( MPI_COMM_WORLD, &rank );
            printf("je suis le processus de rang %d et je mexecute sur %s\n",rank,procname);
        #endif
        error.printErrorStack();
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printErrorStack();
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printErrorStack();
    }
    // catch failure caused by the DataType operations
    catch( DataTypeIException error )
    {
        error.printErrorStack();
    }
};


template <typename TFloeGroup, typename TDynamicsMgr>
void HDF5Manager<TFloeGroup, TDynamicsMgr>::write_boundaries() {
    
    H5File& file( *m_out_file );
    const int   SPACE_DIM = 2;

    Group floe_state_group;
    try {
        floe_state_group = file.openGroup("floe_outlines");
    } catch (...) {
        floe_state_group = file.createGroup("floe_outlines");
    }
    
    for (std::size_t i = 0; i!= m_data_chunk_boundaries.size(); ++i)
    {
        const int   RANK = 3;
        auto& floe_chunk = m_data_chunk_boundaries[i];
        hsize_t     dimsf[RANK] = {m_step_count - m_chunk_step_count, floe_chunk[0].size(), SPACE_DIM};
        hsize_t chunk_dims[RANK] = {m_chunk_step_count, dimsf[1], dimsf[2]};

        DataSet dataset;
        try {
            dataset = floe_state_group.openDataSet(std::to_string(i).c_str());
        } catch (...) {
            FloatType datatype( PredType::NATIVE_DOUBLE );
            datatype.setOrder( H5T_ORDER_LE );
            hsize_t maxdims[RANK] = {H5S_UNLIMITED, dimsf[1], dimsf[2]};
            DataSpace dataspace( RANK, dimsf, maxdims );
            // Modify dataset creation property to enable chunking
            DSetCreatPropList prop;
            prop.setChunk(RANK, chunk_dims);

            dataset = floe_state_group.createDataSet(std::to_string(i).c_str(), datatype, dataspace, prop);
        }

        // Extend the dataset.
        dimsf[0] += m_chunk_step_count;
        dataset.extend(dimsf); 

        DataSpace filespace = dataset.getSpace();
        hsize_t offset[RANK] = {m_step_count - m_chunk_step_count, 0, 0}; // for dataset extension
        filespace.selectHyperslab(H5S_SELECT_SET, chunk_dims, offset);

        // Define memory space.
        DataSpace memspace{RANK, chunk_dims, NULL};

        boost::multi_array<real_type, 3> data(boost::extents[m_chunk_step_count][dimsf[1]][dimsf[2]]);
        for (std::size_t j = 0; j != floe_chunk.size(); ++j)
        {
            for (std::size_t k = 0; k != floe_chunk[j].size(); ++k)
            {
                data[j][k][0] = floe_chunk[j][k][0];
                data[j][k][1] = floe_chunk[j][k][1];
            }
        }

        // Write data to the extended portion of the dataset.
        dataset.write(data.data(), PredType::NATIVE_DOUBLE, memspace, filespace);
    }

    // clearing buffer
    m_data_chunk_boundaries.clear();

};

template <typename TFloeGroup, typename TDynamicsMgr>
void HDF5Manager<TFloeGroup, TDynamicsMgr>::write_states() {
    
    H5File& file( *m_out_file );
    const int   RANK = 3;

    /* saving time */
    DataSet states_dataset;
    const hsize_t nb_floes = m_data_chunk_states[0].size();
    hsize_t     dims[RANK] = {m_step_count - m_chunk_step_count, nb_floes, array_size<saved_state_type>::size };
    const hsize_t     chunk_dims[RANK] = {m_chunk_step_count, dims[1], dims[2]};
    try {
        states_dataset = file.openDataSet("floe_states");
    } catch (...) {
        FloatType datatype( PredType::NATIVE_DOUBLE );
        datatype.setOrder( H5T_ORDER_LE );
        hsize_t maxdims[RANK] = {H5S_UNLIMITED, dims[1], dims[2]}; 
        DataSpace dataspace( RANK, dims, maxdims );
        // Modify dataset creation property to enable chunking
        DSetCreatPropList prop;
        prop.setChunk(RANK, chunk_dims);

        states_dataset = file.createDataSet("floe_states", datatype, dataspace, prop);
    }
    // Extend the dataset.
    dims[0] += chunk_dims[0];
    states_dataset.extend(dims); 

    DataSpace filespace = states_dataset.getSpace();
    hsize_t offset[RANK] = {m_step_count - m_chunk_step_count, 0, 0};
    filespace.selectHyperslab(H5S_SELECT_SET, chunk_dims, offset);
    // Define memory space.
    DataSpace memspace{RANK, chunk_dims, NULL};

    states_dataset.write(m_data_chunk_states.data(), PredType::NATIVE_DOUBLE, memspace, filespace);
};

template <typename TFloeGroup, typename TDynamicsMgr>
void HDF5Manager<TFloeGroup, TDynamicsMgr>::write_time() {
    
    H5File& file( *m_out_file );

    /* saving time */
    DataSet time_dataset;
    hsize_t     dimst[1] = {m_step_count - m_chunk_step_count};
    hsize_t     chunk_dimst[1] = {m_chunk_step_count};
    try {
        time_dataset = file.openDataSet("time");
    } catch (...) {
        FloatType datatype( PredType::NATIVE_DOUBLE );
        datatype.setOrder( H5T_ORDER_LE );
        hsize_t maxdims[1] = {H5S_UNLIMITED}; 
        DataSpace dataspace( 1, dimst, maxdims );
        // Modify dataset creation property to enable chunking
        DSetCreatPropList prop;
        prop.setChunk(1, chunk_dimst);

        time_dataset = file.createDataSet("time", datatype, dataspace, prop);
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

template <typename TFloeGroup, typename TDynamicsMgr>
void HDF5Manager<TFloeGroup, TDynamicsMgr>::write_mass_center() {
    
    H5File& file( *m_out_file );
    const int   RANK = 2;

    /* saving mass center */
    DataSet dataset;
    hsize_t     dimst[RANK] = {m_step_count - m_chunk_step_count, 2};
    hsize_t     chunk_dims[RANK] = {m_chunk_step_count, 2};
    try {
        dataset = file.openDataSet("mass_center");
    } catch (...) {
        FloatType datatype( PredType::NATIVE_DOUBLE );
        datatype.setOrder( H5T_ORDER_LE );
        hsize_t maxdims[RANK] = {H5S_UNLIMITED, 2}; 
        DataSpace dataspace( RANK, dimst, maxdims );
        // Modify dataset creation property to enable chunking
        DSetCreatPropList prop;
        prop.setChunk(RANK, chunk_dims);

        dataset = file.createDataSet("mass_center", datatype, dataspace, prop);
    }
    // Extend the dataset.
    dimst[0] += chunk_dims[0];
    dataset.extend(dimst); 

    DataSpace filespace = dataset.getSpace();
    hsize_t offset[RANK] = {m_step_count - m_chunk_step_count, 0};
    filespace.selectHyperslab(H5S_SELECT_SET, chunk_dims, offset);
    // Define memory space.
    DataSpace memspace{RANK, chunk_dims, NULL};

    dataset.write(m_data_chunk_mass_center.data(), PredType::NATIVE_DOUBLE, memspace, filespace);
};

template <typename TFloeGroup, typename TDynamicsMgr>
void HDF5Manager<TFloeGroup, TDynamicsMgr>::write_OBL_speed() {
    
    H5File& file( *m_out_file );
    const int   RANK = 2;

    /* saving mass center */
    DataSet dataset;
    hsize_t     dimst[RANK] = {m_step_count - m_chunk_step_count, 2};
    hsize_t     chunk_dims[RANK] = {m_chunk_step_count, 2};
    try {
        dataset = file.openDataSet("OBL_speed");
    } catch (...) {
        FloatType datatype( PredType::NATIVE_DOUBLE );
        datatype.setOrder( H5T_ORDER_LE );
        hsize_t maxdims[RANK] = {H5S_UNLIMITED, 2}; 
        DataSpace dataspace( RANK, dimst, maxdims );
        // Modify dataset creation property to enable chunking
        DSetCreatPropList prop;
        prop.setChunk(RANK, chunk_dims);

        dataset = file.createDataSet("OBL_speed", datatype, dataspace, prop);
    }
    // Extend the dataset.
    dimst[0] += chunk_dims[0];
    dataset.extend(dimst); 

    DataSpace filespace = dataset.getSpace();
    hsize_t offset[RANK] = {m_step_count - m_chunk_step_count, 0};
    filespace.selectHyperslab(H5S_SELECT_SET, chunk_dims, offset);
    // Define memory space.
    DataSpace memspace{RANK, chunk_dims, NULL};

    dataset.write(m_data_chunk_OBL_speed.data(), PredType::NATIVE_DOUBLE, memspace, filespace);
};

template <typename TFloeGroup, typename TDynamicsMgr>
void HDF5Manager<TFloeGroup, TDynamicsMgr>::write_kinE() {
    
    H5File& file( *m_out_file );

    /* saving kinE */
    DataSet kinE_dataset;
    hsize_t     dimst[1] = {m_step_count - m_chunk_step_count};
    hsize_t     chunk_dimst[1] = {m_chunk_step_count};
    try {
        kinE_dataset = file.openDataSet("Kinetic Energy");
    } catch (...) {
        FloatType datatype( PredType::NATIVE_DOUBLE );
        // datatype.setOrder( H5T_ORDER_LE );
        hsize_t maxdims[1] = {H5S_UNLIMITED}; 
        DataSpace dataspace( 1, dimst, maxdims );
        // Modify dataset creation property to enable chunking
        DSetCreatPropList prop;
        prop.setChunk(1, chunk_dimst);

        kinE_dataset = file.createDataSet("Kinetic Energy", datatype, dataspace, prop);
    }
    // Extend the dataset.
    dimst[0] += chunk_dimst[0];
    kinE_dataset.extend(dimst); 

    DataSpace filespace = kinE_dataset.getSpace();
    hsize_t offset[1] = {m_step_count - m_chunk_step_count};
    filespace.selectHyperslab(H5S_SELECT_SET, chunk_dimst, offset);
    // Define memory space.
    DataSpace memspace{1, chunk_dimst, NULL};
    kinE_dataset.write(m_data_chunk_kinE, PredType::NATIVE_DOUBLE, memspace, filespace);
};

template <typename TFloeGroup, typename TDynamicsMgr>
double HDF5Manager<TFloeGroup, TDynamicsMgr>::recover_states(
        H5std_string filename, real_type time, floe_group_type& floe_group,
        dynamics_mgr_type& dynamics_manager, bool keep_as_outfile)
{
    
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
    std::vector<real_type> data_time(dims_out[0]);
    for (std::size_t j = 0; j!= dims_out[0]; ++j)
        data_time[j] = 0;
     /*
    * Define the memory dataspace.
    */
    DataSpace time_memspace( 1, dims_out );
    /*
    * Read data from the file
    */
    time_dataset.read( data_time.data(), PredType::NATIVE_DOUBLE, time_memspace, time_dataspace );
    // find time index to read
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
    boost::multi_array<real_type, 2> data_out(boost::extents[dims_out[1]][dims_out[2]]);
    dataset.read( data_out.data(), PredType::NATIVE_DOUBLE, memspace, dataspace );

    std::size_t floe_id = 0;
    for (auto& floe : floe_group.get_floes())
    {
        floe.set_state({
            {data_out[floe_id][0], data_out[floe_id][1]}, data_out[floe_id][2],
            {data_out[floe_id][3], data_out[floe_id][4]}, data_out[floe_id][5],
            {0,0}
        });
        floe.reset_impulse(data_out[floe_id][6]);
        floe_id++;
    }

    {
    // Load OBL speed
    DataSet dataset = file.openDataSet( "OBL_speed" );
    DataSpace dataspace = dataset.getSpace();
    const int rank = 2;
    hsize_t dims_out[rank];
    dataspace.getSimpleExtentDims( dims_out, NULL);
    hsize_t      offset[rank] = {i, 0};  // hyperslab offset in the file
    hsize_t      count[rank] = {1, dims_out[1]};    // size of the hyperslab in the file
    dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );
    hsize_t     dimsm[1] {dims_out[1]};              /* memory space dimensions */
    DataSpace memspace( 1, dimsm );
    std::vector<real_type> data_out(dims_out[1]);
    dataset.read( data_out.data(), PredType::NATIVE_DOUBLE, memspace, dataspace );
    point_type OBL_speed{data_out[0], data_out[1]};
    dynamics_manager.set_OBL_speed(OBL_speed);
    // std::cout << " OBL " << OBL_speed;
    }

    }

    if (keep_as_outfile and i + 1 == dims_out[0]){
        // We keep recover file as output file
        m_step_count = i + 1;
        m_out_file_name = filename;
    }

    return data_time[i];

};


template <
    typename TFloeGroup,
    typename TDynamicsMgr
>
void HDF5Manager<TFloeGroup, TDynamicsMgr>::write_shapes() {
    
    H5File& file( *m_out_file );
    const int   SPACE_DIM = 2;

    /* write shapes */
    try {
        Group floe_shape_group = file.openGroup("floe_shapes");
    } catch (...) {
        /* Create group for floe shapes */
        Group floe_shape_group(file.createGroup("floe_shapes"));

        const int   RANK = 2;
        FloatType datatype( PredType::NATIVE_DOUBLE );
        datatype.setOrder( H5T_ORDER_LE );
        hsize_t     dimsf[2];              // dataset dimensions
        dimsf[1] = SPACE_DIM;
        for (std::size_t i=0; i!=this->nb_considered_floes(); ++i)
        {
            auto& boundary = this->get_floe(i).get_static_floe().geometry().outer();
            dimsf[0] = boundary.size();
            DataSpace dataspace( RANK, dimsf );
            /*
             * Create a nFew dataset within the file using defined dataspace and
             * datatype and default dataset creation properties.
             */
            DataSet dataset = floe_shape_group.createDataSet(H5std_string{std::to_string(i)},datatype, dataspace);
            // auto data = new real_type[dimsf[0] * dimsf[1]];
            boost::multi_array<real_type, 2> data(boost::extents[dimsf[0]][dimsf[1]]);
            for (std::size_t j = 0; j!= dimsf[0]; ++j)
            {
                data[j][0] = boundary[j].x;
                data[j][1] = boundary[j].y;
            }
            /*
             * Write the data to the dataset using default memory space, file
             * space, and transfer properties.
             */
            dataset.write( data.data(), PredType::NATIVE_DOUBLE );
        }
    }

};

template <
    typename TFloeGroup,
    typename TDynamicsMgr
>
void HDF5Manager<TFloeGroup, TDynamicsMgr>::write_window(){
    H5File& file( *m_out_file );

    const int   RANK = 1;
    /* saving ocean window limits */
    hsize_t     dims[RANK] = {4};
    FloatType datatype( PredType::NATIVE_DOUBLE );
    datatype.setOrder( H5T_ORDER_LE );
    DataSpace dataspace( RANK, dims, dims );
    DataSet dataset = file.createDataSet("window", datatype, dataspace);
    auto const& win = m_floe_group->get_initial_window();
    real_type data[4] = {win[0], win[1], win[2], win[3]};
    dataset.write( data, PredType::NATIVE_DOUBLE );
};

template <typename TFloeGroup, typename TDynamicsMgr>
void HDF5Manager<TFloeGroup, TDynamicsMgr>::make_input_file(const dynamics_mgr_type& dynamics_manager){
    try
    {   
        // Prepare manager for writting an input file
        flush();
        save_step(0, dynamics_manager);
        auto saved_out_filename = m_out_file_name;
        /*
         * Turn off the auto-printing when failure occurs so that we can
         * handle the errors appropriately
         */
        Exception::dontPrint();

        // Temporarily change out_file name
        int nb_floe = this->nb_considered_floes();
        int conc = round(m_floe_group->initial_concentration() * 100);
        m_out_file_name = "io/inputs/in_" + std::to_string(nb_floe)
            + "f_" + std::to_string(conc) + "p_" + floe::random::gen_random(5) + ".h5";
        const H5std_string  FILE_NAME( m_out_file_name );

        m_out_file = new H5File( FILE_NAME.c_str(), H5F_ACC_TRUNC );
    
        write_shapes();
        write_window();
        write_states();

        // Close the file after each flush to keep a valid ouput even if program crashes
        delete m_out_file;
        std::cout << m_out_file_name << " written" << std::endl;

        // Reset initial state
        m_out_file_name = saved_out_filename;
        m_chunk_step_count = 0;
        m_step_count = 0;

    }  // end of try block
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printErrorStack();
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printErrorStack();
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printErrorStack();
    }
    // catch failure caused by the DataType operations
    catch( DataTypeIException error )
    {
        error.printErrorStack();
    }
};


template <typename TFloeGroup, typename TDynamicsMgr>
void HDF5Manager<TFloeGroup, TDynamicsMgr>::recover_restrained_floes(const H5std_string filename){

    const H5std_string DATA_SET( "selected_floe_ids" );
    /*
    * Try block to detect exceptions raised by any of the calls inside it
    */
    try{
        /*
         * Turn off the auto-printing when failure occurs so that we can
         * handle the errors appropriately
         */
        Exception::dontPrint();
        H5File* file;
        /*
         * Open the specified file and the specified dataset in the file.
         */
        file = new H5File( filename, H5F_ACC_RDONLY );

        DataSet* dataset = new DataSet(file->openDataSet( DATA_SET ));
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
        int selection_tmp[dims_out[0]];
        dataset->read( selection_tmp, PredType::NATIVE_INT);

        std::vector<std::size_t> selected_floe_ids;
        std::cout << "floe recovered: [";
        for (hsize_t i=0; i<dims_out[0]; ++i) {
            selected_floe_ids.push_back(selection_tmp[i]);
            std::cout << selection_tmp[i] << ", ";
        }
        std::cout << "]" << std::endl;

        this->restrain_floe_ids(selected_floe_ids);
        std::cout << "the selection of floes has been recovered!" << std::endl;
    }
    catch( FileIException error )
    {
        error.printErrorStack();
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printErrorStack();
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printErrorStack();
    }
};

template <typename TFloeGroup, typename TDynamicsMgr>
void HDF5Manager<TFloeGroup, TDynamicsMgr>::write_selected_floe_ids(std::vector<std::size_t> selected_floe_ids){

    const H5std_string FILE_NAME( "io/outputs/selected_floes.h5" );
    const H5std_string SFI("selected_floe_ids");

    /*
     * Try block to detect exceptions raised by any of the calls inside it
     */
    try{
        /*
         * Turn off the auto-printing when failure occurs so that we can
         * handle the errors appropriately
         */
        Exception::dontPrint();
        /*
         * Create a file.
         */
        H5File* file;
        file = new H5File( FILE_NAME, H5F_ACC_TRUNC );

        /* write list of selected floes */
        const hsize_t dim[1] = {selected_floe_ids.size()};
        DataSpace space( 1, dim );

        DataSet* data_floes = new DataSet(file->createDataSet( SFI, PredType::NATIVE_INT, space ));

        int val[selected_floe_ids.size()];
        for (size_t i=0; i<selected_floe_ids.size(); ++i) {
            val[i] = selected_floe_ids[i];
        }

        data_floes->write(val, PredType::NATIVE_INT);

        delete data_floes;
        delete file;
    } // end of try block
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printErrorStack();
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printErrorStack();
    }
};

}} // namespace floe::io


// #endif // FLOE_IO_HDF5_MANAGER_DEF_HPP
