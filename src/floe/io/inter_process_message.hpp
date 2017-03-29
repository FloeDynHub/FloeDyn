/*!
 * \file floe/io/inter_process_message.hpp
 * \brief
 * \author Quentin Jouet
 */

#ifndef FLOE_IO_INTER_PROCESS_MESSAGE_HPP
#define FLOE_IO_INTER_PROCESS_MESSAGE_HPP

#include <map>
#include <list>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/array.hpp>


namespace floe { namespace io
{

enum JobTag {
    collision_job,
    time_step_job,
    move_job,
    interpene_job,
    test_job,
    termination_signal
};

template <
    typename T
>
class InterProcessMessage
{

public:

    // Type traits
    using value_type = T;
    using id_list_type = std::vector<std::size_t>;

    //! Default constructor.
    InterProcessMessage() {}
    //! Constructor with id
    InterProcessMessage(int id) : m_id{id} {}
    //! Constructor from request (for response)
    InterProcessMessage(InterProcessMessage const& request) :
        m_id{request.id()}, m_tag{request.tag()} {}

    //! Accessors
    inline std::map<int, std::array<value_type, 7>> const&  states() const { return m_states; }
    inline id_list_type const&  floe_ids() const { return m_floe_ids; }
    inline void set_floe_ids(id_list_type const& floe_ids) { m_floe_ids = floe_ids; }
    inline JobTag tag() const { return m_tag; }
    inline void set_tag(JobTag tag) { m_tag = tag; }
    inline int id() const { return m_id; }
    inline void store_time_step(value_type delta_t) { m_delta_t = delta_t; }
    inline value_type time_step() const { return m_delta_t; }
    inline void store_time(value_type t) { m_time = t; }
    inline value_type time() const { return m_time; }
    inline int nb_LCP_solved() const { return m_nb_LCP_solved; }
    inline void nb_LCP_solved(int n) { m_nb_LCP_solved = n; }
    inline bool interpenetration() const { return m_interpenetration; }
    inline void interpenetration(bool b) { m_interpenetration = b; }
    inline int mpi_source() const { return m_mpi_source; }
    inline void mpi_source(int n) { m_mpi_source = n; }

    //! Setter
    template<typename TFloeGroup, typename IdsIterable>
    void store_states(TFloeGroup const& floe_group, IdsIterable const& floe_ids){
        m_floe_ids = floe_ids;
        for (std::size_t id : floe_ids){
            auto const& floe = floe_group.get_floes()(id);
            auto const& state = floe.state();
            m_states[id] = {{
                state.pos.x, state.pos.y, state.theta,
                state.speed.x, state.speed.y, state.rot,
                floe.total_received_impulse()
            }};
        }
    }

    template<typename TFloeGroup, typename IdsIterable>
    void store_states_light(TFloeGroup const& floe_group, IdsIterable const& floe_ids, int mpi_dest){
        m_floe_ids = floe_ids;
        for (int id : floe_ids){
            if (floe_group.states_origin()[id] != mpi_dest){
                auto const& floe = floe_group.get_floes()(id);
                auto const& state = floe.state();
                m_states[id] = {{
                    state.pos.x, state.pos.y, state.theta,
                    state.speed.x, state.speed.y, state.rot,
                    floe.total_received_impulse()
                }};
            }
        }
    }

    // This method lets cereal know which data members to serialize
    template<class Archive>
    void serialize(Archive & archive)
    {
    archive( m_id, m_tag, m_floe_ids, m_states, m_delta_t,
        m_time, m_nb_LCP_solved, m_interpenetration ); // serialize things by passing them to the archive
    }

private:
    int m_id;
    JobTag m_tag;
    id_list_type m_floe_ids;
    std::map<int, std::array<value_type, 7>> m_states;
    value_type m_delta_t = 0;
    value_type m_time;
    int m_nb_LCP_solved = 0;
    bool m_interpenetration = false;
    int m_mpi_source = -1;
};

}} // namespace floe::io


#endif // FLOE_IO_INTER_PROCESS_MESSAGE_HPP

