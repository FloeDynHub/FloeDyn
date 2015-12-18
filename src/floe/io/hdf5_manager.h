/*!
 * \file io/hdf5_manager.h
 * \brief HDF5 manager for io
 * \author Quentin Jouet
 */

#ifndef FLOE_IO_HDF5_MANAGER_HPP
#define FLOE_IO_HDF5_MANAGER_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include "floe/floes/floe_group.hpp"

#include "H5Cpp.h"
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif


namespace floe { namespace io
{

/*! HDF5Manager
 *
 * Handles :
 *   Floe outlines, floe states, ocean state, and time output
 *   Floe states input
 *
 * Internally saves a fixed number of states before writing it to the file.
 *
 */

// Generate random string
std::string gen_random(const int len);


// Get array size at compile time
template<typename>
struct array_size;

template<typename T, size_t N>
struct array_size<std::array<T,N> > {
    static size_t const size = N;
};


template <typename TFloeGroup, typename TDynamicsMgr>
class HDF5Manager
{

public:
    template<typename T>
    using vector = std::vector<T>;
    using floe_group_type = TFloeGroup;
    using dynamics_mgr_type = TDynamicsMgr;
    using value_type = typename TFloeGroup::value_type;
    using point_type = typename TFloeGroup::point_type; 
    using saved_state_type = std::array<value_type, 9>;

    //! Default constructor.
    HDF5Manager(floe_group_type const& floe_group);
    //! Destructor
    ~HDF5Manager();

    //! Save the current simulation state for output
    void save_step(value_type time, const dynamics_mgr_type&);
    //! Recover simulation state from file
    double recover_states(H5std_string filename, value_type time, floe_group_type&, dynamics_mgr_type&);
    //! Recover simulation state from file
    inline void set_floe_group(floe_group_type const& floe_group) { m_floe_group = &floe_group; };

private:

    std::string m_out_file_name; //!< output file name
    H5File* m_out_file; //!< output file
    hsize_t m_step_count; //!< Total nb of outputted simulation states
    hsize_t m_chunk_step_count; //! Nb of temporarily saved steps (to flush in out file)
    const hsize_t m_flush_max_step; //! Max nb of temporarily saved steps (chunk size)
    floe_group_type const* m_floe_group;

    vector<vector<vector<vector<value_type>>>> m_data_chunk_boundaries; //!< Temp saved floe boundaries
    vector<vector<saved_state_type>> m_data_chunk_states; //!< Temp saved floe states
    value_type* m_data_chunk_time; //!< Temp saved times
    vector<std::array<value_type, 2>> m_data_chunk_mass_center; //!< Temp saved floe group mass centers
    vector<std::array<value_type, 2>> m_data_chunk_OBL_speed; //!< Temp saved ocean datas

    //! Flush temporarily saved data
    void write_chunk();

    //! out floe shapes (boundary in relative frame)
    void write_shapes();

    //! Partial writings :
    void write_boundaries();
    void write_states();
    void write_time();
    void write_mass_center();
    void write_OBL_speed();

};


}} // namespace floe::io


#endif // FLOE_IO_HDF5_MANAGER_HPP
