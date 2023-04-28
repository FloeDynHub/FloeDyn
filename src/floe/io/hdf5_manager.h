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
#include <math.h>
#include <memory>
#include "boost/multi_array.hpp"
#include "floe/floes/floe_group.hpp"

#include "H5Cpp.h"
// #ifndef H5_NO_NAMESPACE
// using namespace H5;
// #endif


namespace floe { namespace io
{
    using namespace H5;

/*! HDF5Manager
 *
 * Handles :
 *   Floe outlines, floe states, ocean state, and time output
 *   Floe states input
 *
 * Internally saves a fixed number of states before writing it to the file.
 *
 */

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
    using real_type = typename TFloeGroup::real_type;
    using point_type = typename TFloeGroup::point_type; 
    using saved_state_type = std::array<real_type, 10>;

    //! Default constructor.
    HDF5Manager(floe_group_type const& floe_group);
    //! Destructor
    ~HDF5Manager();

    //! make input file from floe_group
    void make_input_file(const dynamics_mgr_type& dynamics_manager);
    //! Calls save_step() if this time needs to be saved
    void save_step_if_needed(real_type time, const dynamics_mgr_type&);
    //! Save the current simulation state for output
    void save_step(real_type time, const dynamics_mgr_type&);
    //! Flush temporarily saved data
    void flush();
    //! Recover simulation state from file
    double recover_states(H5std_string filename, real_type time, floe_group_type&,
                          dynamics_mgr_type&, bool keep_as_outfile);
    void set_floe_group(floe_group_type const& floe_group) {
        m_floe_group = &floe_group; 
        m_data_chunk_states.resize(boost::extents[m_flush_max_step][this->nb_considered_floes()][10]);
    };
    //! Do not consider all floes in floe group, /!\ only do this at the begining (resizes out states dataset)
    void restrain_floe_ids(std::vector<std::size_t> id_list) {
        m_floe_ids = id_list; 
        m_data_chunk_states.resize(boost::extents[m_flush_max_step][m_floe_ids.size()][10]);
    };
    inline bool is_restrained() const { return m_floe_ids.size(); }
    inline std::string const& out_file_name() const { return m_out_file_name; }
    inline void set_out_file_name(std::string file_name) { m_out_file_name = file_name; }
    inline real_type get_next_out_limit() const {return m_next_out_limit;}
    inline void set_out_step(real_type out_step, real_type time) {
        m_out_step = out_step;
        if (time > 0) m_next_out_limit = (std::floor(time / m_out_step) + 1) * m_out_step;
    }
    inline void auto_step_count(real_type time){ 
        this->m_step_count = (int)time/this->m_out_step; 
        std::cout << "It is getting back the storage at the counter step: " << this->m_step_count << "\n";
    }

    void recover_restrained_floes(const H5std_string filename);

    void write_selected_floe_ids(std::vector<std::size_t> selected_floe_ids);



private:

    std::string m_out_file_name; //!< output file name
    H5File* m_out_file; //!< output file
    Group* m_shapes_group;
    hsize_t m_step_count; //!< Total nb of outputted simulation states
    hsize_t m_chunk_step_count; //!< Nb of temporarily saved steps (to flush in out file)
    const hsize_t m_flush_max_step; //!< Max nb of temporarily saved steps (chunk size)
    floe_group_type const* m_floe_group; //!< floe group pointer
    std::vector<std::size_t> m_floe_ids; //!< Restrained id list to consider for output

    vector<vector<vector<vector<real_type>>>> m_data_chunk_boundaries; //!< Temp saved floe boundaries
    boost::multi_array<real_type, 3> m_data_chunk_states; //!< Temp saved floe states
    real_type* m_data_chunk_time; //!< Temp saved times
    boost::multi_array<real_type, 2> m_data_chunk_mass_center; //!< Temp saved floe group mass centers
    boost::multi_array<real_type, 2> m_data_chunk_OBL_speed; //!< Temp saved ocean datas
    real_type* m_data_chunk_kinE; //!< Temp saved Kinetic Energy

    // output
    real_type m_out_step; //!< Time step between simulation state outputs
    real_type m_next_out_limit; //!< Next time limit for state ouput
    //! Number of floe shapes written to file (fracture creates new ones)
    hsize_t m_nb_floe_shapes_written;

    //! out floe shapes (boundary in relative frame)
    void write_shapes();
    //! Partial writings :
    void write_boundaries();
    void write_states();
    void write_time();
    void write_mass_center();
    void write_OBL_speed();
    void write_window();
    void write_kinE();

    inline std::size_t nb_considered_floes() const { return m_floe_ids.size() ? m_floe_ids.size() : m_floe_group->get_floes().size(); }
    inline typename floe_group_type::floe_type const& get_floe(std::size_t id) const {
        return m_floe_ids.size() ? m_floe_group->get_floes()[m_floe_ids[id]] : m_floe_group->get_floes()[id];
    }

    inline void update_next_out_limit() { m_next_out_limit += m_out_step; }
    inline bool need_step_output(real_type time) { return (m_out_step && time >= m_next_out_limit); }

};


}} // namespace floe::io


#endif // FLOE_IO_HDF5_MANAGER_HPP
