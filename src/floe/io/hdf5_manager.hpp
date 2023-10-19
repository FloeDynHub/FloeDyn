/*!
 * \file floes/hdf5_manager.hpp
 * \brief HDF5 manager for io
 * \author Quentin Jouet
 */

#ifndef FLOE_IO_HDF5_MANAGER_DEF_HPP
#define FLOE_IO_HDF5_MANAGER_DEF_HPP

#include "floe/io/hdf5_manager.h"
#include "floe/utils/random.hpp"

#define WHEREAMI std::cout << std::endl << "no crash until line " << __LINE__ << " in the file " __FILE__ << std::endl;


namespace floe { namespace io
{

//! Default constructor.
template <typename TFloeGroup, typename TDynamicsMgr>
HDF5Manager<TFloeGroup, TDynamicsMgr>::HDF5Manager(floe_group_type const& floe_group, bool export_mesh) :
    m_out_file_name{"io/outputs/0_test.h5"},
    // m_out_file_name{"io/outputs/out_mesh_" + floe::random::gen_random(5) + ".h5"},
    m_out_file{nullptr}, m_step_count{0}, m_chunk_step_count{0},
    m_flush_max_step{2}, // min val = 2
    m_floe_group{&floe_group},
    m_data_chunk_states(boost::extents[0][0][0]),
    m_data_chunk_elem_data(boost::extents[0][0][0]),
    m_data_chunk_node_data(boost::extents[0][0][0]),
    m_data_chunk_time{new real_type[m_flush_max_step]},
    m_data_chunk_mass_center(boost::extents[m_flush_max_step][2]),
    m_data_chunk_OBL_speed(boost::extents[m_flush_max_step][2]),
    m_data_chunk_kinE{new real_type[m_flush_max_step]},
    m_out_step{0}, m_next_out_limit{0}, m_nb_floe_shapes_written{0}, m_nb_floe_meshes_coord_written{0}, m_nb_floe_meshes_connect_written{0}, m_shapes_group{nullptr}, m_meshes_coord_group{nullptr}, m_meshes_connect_group{nullptr},
    m_max_elem{0},
    m_max_nodes{0},
    m_export_mesh{export_mesh}
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
    m_max_elem = m_floe_group->get_max_elem();
    m_max_nodes = m_floe_group->get_max_nodes();
            
    // Handle fracture
    if (m_data_chunk_states.size() > 0 && m_data_chunk_states[0].size() != this->nb_considered_floes()) {
        flush();
        m_chunk_step_count = 0;
        m_data_chunk_states.resize(boost::extents[m_flush_max_step][this->nb_considered_floes()][array_size<saved_state_type>::size]);
        m_data_chunk_elem_data.resize(boost::extents[m_flush_max_step][this->nb_considered_floes()][m_max_elem]);
        m_data_chunk_node_data.resize(boost::extents[m_flush_max_step][this->nb_considered_floes()][m_max_nodes]);
        write_shapes();
        if (m_export_mesh)
        {
            write_meshes_coord();
            write_meshes_connect();
        }
    }
    // save states
    if (m_data_chunk_states.size() == 0) m_data_chunk_states.resize(boost::extents[m_flush_max_step][this->nb_considered_floes()][array_size<saved_state_type>::size]);
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
            floe.state().pos.y,
            (floe.state().is_active()) ? 1. : 0.,
            floe.static_floe().thickness()
            
        }){
            m_data_chunk_states[m_chunk_step_count][id][k++] = val;
        }
    }
    if (m_export_mesh)
    {
        // save elem data
        if (m_data_chunk_elem_data.size() == 0) m_data_chunk_elem_data.resize(boost::extents[m_flush_max_step][this->nb_considered_floes()][m_max_elem]); 
        for(std::size_t iFloe = 0; iFloe < this->nb_considered_floes(); ++iFloe)
        {
            auto const& floe = this->get_floe(iFloe);
            std::vector<std::vector<real_type>> femSolStress = floe.get_fem_stress();
            for (std::size_t iElem = 0 ; iElem < floe.mesh().get_n_cells() ; ++iElem)
            {
                if (femSolStress.size() == floe.mesh().get_n_cells())
                {
                    // std::cout << "Writing stress to hdf5 output file" << std::endl ;
                    m_data_chunk_elem_data[m_chunk_step_count][iFloe][iElem] = (real_type)sqrt(pow(femSolStress[iElem][0], 2) + pow(femSolStress[iElem][1], 2) - femSolStress[iElem][0]*femSolStress[iElem][1] + 3*pow(femSolStress[iElem][2], 2));
                }
                else
                {
                    // std::cout << "Nope." << std::endl ;
                    m_data_chunk_elem_data[m_chunk_step_count][iFloe][iElem] = (real_type)iElem;
                    // m_data_chunk_elem_data[m_chunk_step_count][iFloe][iElem] = floe.total_received_impulse();
                }
            }
        }
        // saving nodal data 
        if (m_data_chunk_node_data.size() == 0) m_data_chunk_node_data.resize(boost::extents[m_flush_max_step][this->nb_considered_floes()][m_max_nodes]); 
        for(std::size_t iFloe = 0; iFloe < this->nb_considered_floes(); ++iFloe)
        {
            auto const& floe = this->get_floe(iFloe);
            std::vector<real_type> femSol = floe.get_fem_solution();
            if (femSol.size() != floe.mesh().get_n_nodes())
            {
                std::cout << "Incoherent size. Size of femSol : " << femSol.size() << " instead of " << floe.mesh().get_n_nodes() << std::endl;
            }

            for (std::size_t iNode = 0 ; iNode < floe.mesh().get_n_nodes() ; ++iNode)
            {
                // m_data_chunk_node_data[m_chunk_step_count][iFloe][iNode] = (real_type)iNode;
                if (femSol.size() == floe.mesh().get_n_nodes())
                {
                    // std::cout << "I got something interesting to write !! " << std::endl;
                    m_data_chunk_node_data[m_chunk_step_count][iFloe][iNode] = (real_type)femSol[iNode];
                }
                else 
                    m_data_chunk_node_data[m_chunk_step_count][iFloe][iNode] = (real_type)iNode;

            }
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
        catch (...) {
            m_shapes_group = new Group( m_out_file->createGroup("floe_shapes") );
            write_shapes();
        }
        if (m_export_mesh)
        {
            try { m_out_file->openGroup("floe_meshes_coord"); }
            catch (...) {
                m_meshes_coord_group = new Group( m_out_file->createGroup("floe_meshes_coord") );
                write_meshes_coord();
            }
            try { m_out_file->openGroup("floe_meshes_connect"); }
            catch (...) {
                m_meshes_connect_group = new Group( m_out_file->createGroup("floe_meshes_connect") );
                write_meshes_connect();
            }
        }

        try { m_out_file->openDataSet("window"); }
        catch (...) { write_window(); }

        // write_boundaries();
        write_states();
        if (m_export_mesh)
        {
            write_elem_data();
            write_node_data();
        }
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
    hsize_t     dims[RANK] = {m_step_count - m_chunk_step_count, nb_floes, array_size<saved_state_type>::size};
    const hsize_t     chunk_dims[RANK] = {m_chunk_step_count, dims[1], dims[2]};
    try {
        states_dataset = file.openDataSet("floe_states");
    } catch (...) {
        FloatType datatype( PredType::NATIVE_DOUBLE );
        datatype.setOrder( H5T_ORDER_LE );
        hsize_t maxdims[RANK] = {H5S_UNLIMITED, H5S_UNLIMITED, dims[2]}; 
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
void HDF5Manager<TFloeGroup, TDynamicsMgr>::write_elem_data() {
    
    H5File& file( *m_out_file );
    const int   RANK = 3;

    /* saving time */
    DataSet elem_data_dataset;
    const hsize_t nb_floes = m_data_chunk_elem_data[0].size();

    m_max_elem = m_floe_group->get_max_elem();

    hsize_t     dims[RANK] = {m_step_count - m_chunk_step_count, nb_floes, m_max_elem};
    const hsize_t     chunk_dims[RANK] = {m_chunk_step_count, dims[1], dims[2]};
    try {
        elem_data_dataset = file.openDataSet("floe_elem_data");
    } catch (...) {
        FloatType datatype( PredType::NATIVE_DOUBLE );
        datatype.setOrder( H5T_ORDER_LE );
        hsize_t maxdims[RANK] = {H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED}; 
        DataSpace dataspace( RANK, dims, maxdims );
        // Modify dataset creation property to enable chunking
        DSetCreatPropList prop;
        prop.setChunk(RANK, chunk_dims);

        elem_data_dataset = file.createDataSet("floe_elem_data", datatype, dataspace, prop);
    }
    // Extend the dataset.
    dims[0] += chunk_dims[0];
    elem_data_dataset.extend(dims); 
    DataSpace filespace = elem_data_dataset.getSpace();
    hsize_t offset[RANK] = {m_step_count - m_chunk_step_count, 0, 0};
    filespace.selectHyperslab(H5S_SELECT_SET, chunk_dims, offset);
    // Define memory space.
    DataSpace memspace{RANK, chunk_dims, NULL};

    elem_data_dataset.write(m_data_chunk_elem_data.data(), PredType::NATIVE_DOUBLE, memspace, filespace);
};


template <typename TFloeGroup, typename TDynamicsMgr>
void HDF5Manager<TFloeGroup, TDynamicsMgr>::write_node_data() {
    
    H5File& file( *m_out_file );
    // const int   RANK = 4;
    const int   RANK = 3;
    
    /* saving time */
    DataSet node_data_dataset;
    const hsize_t nb_floes = m_data_chunk_elem_data[0].size();
    // const hsize_t nb_data = array_size<saved_elem_data_type>::size;

    m_max_nodes = m_floe_group->get_max_nodes();

    hsize_t     dims[RANK] = {m_step_count - m_chunk_step_count, nb_floes, m_max_nodes};
    const hsize_t     chunk_dims[RANK] = {m_chunk_step_count, dims[1], dims[2]};
    try {
        node_data_dataset = file.openDataSet("floe_node_data");
    } catch (...) {
        FloatType datatype( PredType::NATIVE_DOUBLE );
        datatype.setOrder( H5T_ORDER_LE );
        hsize_t maxdims[RANK] = {H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED}; 
        DataSpace dataspace( RANK, dims, maxdims );
        // Modify dataset creation property to enable chunking
        DSetCreatPropList prop;
        prop.setChunk(RANK, chunk_dims);
        node_data_dataset = file.createDataSet("floe_node_data", datatype, dataspace, prop);
    }
    // Extend the dataset.
    dims[0] += chunk_dims[0];
    node_data_dataset.extend(dims); 

    DataSpace filespace = node_data_dataset.getSpace();
    hsize_t offset[RANK] = {m_step_count - m_chunk_step_count, 0, 0};
    filespace.selectHyperslab(H5S_SELECT_SET, chunk_dims, offset);
    // Define memory space.
    DataSpace memspace{RANK, chunk_dims, NULL};

    node_data_dataset.write(m_data_chunk_node_data.data(), PredType::NATIVE_DOUBLE, memspace, filespace);
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

    floe_group.get_floes().filter_off(); // crack version
    floe_group.get_floes().resize(dims_out[1]); // Resize for crack version
    floe_group.get_floe_group_h().m_list_floe_h.resize(dims_out[1]);
    for (std::size_t floe_id = 0; floe_id < floe_group.get_floes().size(); floe_id++)
    {
        auto& floe = floe_group.get_floes()[floe_id];
        // auto& new_state = floe.state();
        bool active = (data_out[floe_id][9] == 1.) ? true : false; // crack version
        floe.state().set_active(active); // crack version
        if (active) {
            floe.set_state({
                {data_out[floe_id][0], data_out[floe_id][1]}, data_out[floe_id][2],
                {data_out[floe_id][3], data_out[floe_id][4]}, data_out[floe_id][5],
                {0,0}
            });
            floe.reset_impulse(data_out[floe_id][6]);
            floe.static_floe().set_thickness(data_out[floe_id][10]);
        }
    }
    floe_group.update_list_ids_active(); // crack version
    
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

    const int   RANK = 2;
    FloatType datatype( PredType::NATIVE_DOUBLE );
    datatype.setOrder( H5T_ORDER_LE );
    hsize_t     dimsf[2];              // dataset dimensions
    dimsf[1] = SPACE_DIM;
    for (std::size_t i=m_nb_floe_shapes_written; i!=this->nb_considered_floes(); ++i)
    {
        auto& boundary = this->get_floe(i).get_static_floe().geometry().outer();
        dimsf[0] = boundary.size();
        DataSpace dataspace( RANK, dimsf );
        /*
            * Create a nFew dataset within the file using defined dataspace and
            * datatype and default dataset creation properties.
            */
        DataSet dataset = m_shapes_group->createDataSet(H5std_string{std::to_string(i)},datatype, dataspace);
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
        // add atribute for oceanic skin drag
        DataSpace att_space(H5S_SCALAR);
        auto val = this->get_floe(i).get_static_floe().C_w();
        Attribute att = dataset.createAttribute("C_w", datatype, att_space );
        att.write( datatype, &val );
    }

    m_nb_floe_shapes_written = this->nb_considered_floes();

};


template <
    typename TFloeGroup,
    typename TDynamicsMgr
>
void HDF5Manager<TFloeGroup, TDynamicsMgr>::write_meshes_coord() {
    
    H5File& file( *m_out_file );
    const int   SPACE_DIM = 2;

    const int   RANK = 2;
    FloatType datatype( PredType::NATIVE_DOUBLE );
    datatype.setOrder( H5T_ORDER_LE );
    hsize_t     dimsf[2];              // dataset dimensions
    dimsf[1] = SPACE_DIM;
    for (std::size_t iFloe=m_nb_floe_meshes_coord_written; iFloe != this->nb_considered_floes(); ++iFloe)
    {
        // exporting the coordinate table 
        boost::geometry::model::multi_point<point_type> coord = this->get_floe(iFloe).get_mesh().points();

        dimsf[0] = coord.size();
        DataSpace dataspace( RANK, dimsf );
        DataSet dataset = m_meshes_coord_group->createDataSet(H5std_string{std::to_string(iFloe)},datatype, dataspace);
        boost::multi_array<real_type, 2> data(boost::extents[dimsf[0]][dimsf[1]]);
        for (std::size_t iPoint = 0; iPoint < dimsf[0]; ++iPoint)
        {
            // coordinates are first put in the floe coordinate system. (< -pos.x -pos.y > translation, and -theta rotation)
            real_type x = coord[iPoint][0] - this->get_floe(iFloe).state().pos.x;
            real_type y = coord[iPoint][1] - this->get_floe(iFloe).state().pos.y;
            real_type theta = this->get_floe(iFloe).state().theta;

            data[iPoint][0] = x*cos(theta) + y*sin(theta);
            data[iPoint][1] = -x*sin(theta) + y*cos(theta);
        }
        /*
            * Write the data to the dataset using default memory space, file
            * space, and transfer properties.
            */
        dataset.write( data.data(), PredType::NATIVE_DOUBLE );
        // add attribute for floe id
        DataSpace att_space(H5S_SCALAR);
        auto val = (double)iFloe ;
        Attribute att = dataset.createAttribute("index", datatype, att_space );
        att.write( datatype, &val );
    }
    m_nb_floe_meshes_coord_written = this->nb_considered_floes();
};


template <
    typename TFloeGroup,
    typename TDynamicsMgr
>
void HDF5Manager<TFloeGroup, TDynamicsMgr>::write_meshes_connect() {
    
    H5File& file( *m_out_file );
    // const int   SPACE_DIM = 2;
    const int   RANK = 2;
    FloatType datatype( PredType::NATIVE_DOUBLE );
    datatype.setOrder( H5T_ORDER_LE );
    hsize_t     dimsf[2];              // dataset dimensions
    dimsf[1] = 3; // we're supposed to have only linear triangles. Would'nt it be nice, to, add, a check here... 
    // or even better, add something like this->get_floe(iFloe).get_mesh().getn_max_nodes() 
    for (std::size_t iFloe=m_nb_floe_meshes_connect_written; iFloe != this->nb_considered_floes(); ++iFloe)
    {
        // exporting the connectivity table 
        std::vector<std::array<std::size_t,3>> connect =this->get_floe(iFloe).get_mesh().connectivity(); 
        dimsf[0] = connect.size();
        DataSpace dataspace( RANK, dimsf );
        DataSet dataset = m_meshes_connect_group->createDataSet(H5std_string{std::to_string(iFloe)},datatype, dataspace);
        boost::multi_array<int, 2> data(boost::extents[dimsf[0]][dimsf[1]]);
        for (std::size_t iElem = 0; iElem < dimsf[0]; ++iElem)
        {
            data[iElem][0] = (int)connect[iElem][0];
            data[iElem][1] = (int)connect[iElem][1];
            data[iElem][2] = (int)connect[iElem][2];
        }
        dataset.write( data.data(), PredType::NATIVE_INT );
        // add attribute for floe id
        DataSpace att_space(H5S_SCALAR);
        auto val = (double)iFloe ;
        Attribute att = dataset.createAttribute("index", datatype, att_space );
        att.write( datatype, &val );
    }
    m_nb_floe_meshes_connect_written = this->nb_considered_floes();
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
        m_shapes_group = new Group( m_out_file->createGroup("floe_shapes") );
    
        write_shapes();
        write_window();
        write_states();
        if (m_export_mesh)
        {
            write_elem_data();
            write_node_data();
        }

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


#endif // FLOE_IO_HDF5_MANAGER_DEF_HPP
