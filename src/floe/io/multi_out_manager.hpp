/*!
 * \file floes/multi_out_manager.hpp
 * \brief Wrapper for multiple HDF5 manager for io
 * \author Quentin Jouet
 */

#ifndef FLOE_IO_MULTI_OUT_MANAGER_HPP
#define FLOE_IO_MULTI_OUT_MANAGER_HPP

#include <random>
#include <vector>
#include <type_traits>


namespace floe { namespace io
{

template <typename TOutManager>
class MultiOutManager
{

public:
    using floe_group_type = typename TOutManager::floe_group_type;
    using dynamics_mgr_type = typename TOutManager::dynamics_mgr_type;
    using real_type = typename TOutManager::real_type;
    //! Default constructor.
    MultiOutManager(floe_group_type const& floe_group) :
        m_out_managers{floe_group, floe_group},
        m_floe_group{&floe_group}
    {
        TOutManager& second_out_mgr = m_out_managers[1];
        second_out_mgr.set_out_step(10, 0);
        second_out_mgr.set_out_file_name("/Users/matthiasrabatel/Travail/outputs_mycode/out_partial.h5"); // todo : less hardcoded if needed...
    }

    void restrain_second_mgr(){
        // Configure second out manager to manage only a subgroup of floes
        // little out_step + big number of floes = too big out file
        // -> first manager outs all floes rarely, second outs few floes often
        auto const& floe_group = *m_floe_group;
        TOutManager& second_out_mgr = m_out_managers[1];
        // get floe_group window
        auto win = floe_group.get_initial_window();
        real_type x_margin = (win[1] - win[0]) / 5;
        real_type y_margin = (win[3] - win[2]) / 5;
        typename std::remove_reference<decltype(win)>::type subwin{{
            win[0] + x_margin,
            win[1] - x_margin,
            win[2] + y_margin,
            win[3] - y_margin
        }};
        // subgroup of floes in subwin
        std::vector<std::size_t> interesting_floe_ids;
        std::size_t id = 0;
        for (auto const& floe : floe_group.get_floes()){
            if (subwin[0] < floe.state().pos.x
                and floe.state().pos.x < subwin[1]
                and subwin[2] < floe.state().pos.y
                and floe.state().pos.y < subwin[3])
            {
                interesting_floe_ids.push_back(id);
            }
            id++;
        }
        std::srand(0); // be sure to get the same result.
        std::random_shuffle(
            interesting_floe_ids.begin(),
            interesting_floe_ids.end()
        );
        // selecting floes for output
        std::size_t nb_floe_select = 1000;
        std::vector<std::size_t> selected_floe_ids(
            interesting_floe_ids.begin(),
            interesting_floe_ids.begin() + nb_floe_select
        );
        second_out_mgr.restrain_floe_ids(selected_floe_ids);
        // second_out_mgr.set_out_file_name("");
    }

    //! make input file from floe_group
    void make_input_file(const dynamics_mgr_type& dynamics_manager){
        this->m_out_managers[0].make_input_file(dynamics_manager);
    };
    //! Calls save_step() if this time needs to be saved
    void save_step_if_needed(real_type time, const dynamics_mgr_type& dyn_mgr){
        if (!this->m_out_managers[1].is_restrained()){
            this->restrain_second_mgr();
        }
        for (auto& mgr : this->m_out_managers) mgr.save_step_if_needed(time, dyn_mgr);
    }
    //! Flush temporarily saved data
    void flush(){
        for (auto& mgr : this->m_out_managers) mgr.flush();
    }
    //! Recover simulation state from file
    double recover_states(
        H5std_string filename, real_type time, floe_group_type& floe_group,
        dynamics_mgr_type& dyn_mgr, bool keep_as_outfile)
    {
        real_type recover_time = m_out_managers[0].recover_states(
            filename, time, floe_group, dyn_mgr, keep_as_outfile
        );
        m_out_managers[1].auto_step_count(recover_time);
        return recover_time;
    }
    void set_floe_group(floe_group_type const& floe_group) {
        for (auto& mgr : this->m_out_managers) mgr.set_floe_group(floe_group);
    };
    inline std::string const& out_file_name() const { 
        return m_out_managers[0].out_file_name();
    }
    inline void set_out_step(real_type out_step, real_type time) { 
        this->m_out_managers[0].set_out_step(out_step, time);
    }

private:
    std::vector<TOutManager> m_out_managers;
    floe_group_type const* m_floe_group; //!< floe group pointer

};

}} // namespace floe::io


#endif // FLOE_IO_MULTI_OUT_MANAGER_HPP
