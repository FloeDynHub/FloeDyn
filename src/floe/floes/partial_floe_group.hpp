/*!
 * \file floe/floes/partial_floe_group.hpp
 * \brief Partial Floe Configuration class (subset of floes)
 * \author Quentin Jouet
 */

#ifndef FLOES_PARTIAL_FLOE_GROUP_HPP
#define FLOES_PARTIAL_FLOE_GROUP_HPP

#include "floe/floes/floe_group.hpp"
#include "floe/arithmetic/filtered_container.hpp"
#include "floe/io/inter_process_message.hpp"
 
namespace floe { namespace floes
{

/*! FloeGroup
 *
 * It represents a subset of the set of Floes.
 *
 */


template <
    typename TFloe,
    typename TFloeList = FilteredVector<TFloe>
>
class PartialFloeGroup : public FloeGroup<TFloe, TFloeList>
{

public:
    using base_class = FloeGroup<TFloe, TFloeList>;
    using floe_type = TFloe;
    using real_type = typename floe_type::real_type;
    using message_type = io::InterProcessMessage<real_type>;

    void update_partial_list(std::vector<std::size_t> floe_id_list){
        base_class::get_floes().update_ids(floe_id_list);
    }
    virtual int absolute_id(int id) const override {
        return base_class::get_floes().absolute_id(id);
    }
    inline std::vector<int> const& states_origin() const { return m_states_origin; };
    void update_floe_states(message_type const& msg, bool update=true); // override;
    virtual void post_load_floe() override { m_states_origin.clear(); m_states_origin.resize(this->get_floes().size(), 0); }
    virtual void recover_previous_step_states() override { base_class::recover_previous_step_states(); this->post_load_floe(); };
private:
    std::vector<int> m_states_origin;
};


template <typename TFloe, typename TFloeList>
void
PartialFloeGroup<TFloe, TFloeList>::update_floe_states(message_type const& msg, bool update)
{
    for (auto const& iter : msg.states()){
        auto const& s = iter.second;
        auto& floe = this->m_list_floe(iter.first);
        if (update){
            floe.set_state({{s[0], s[1]}, s[2], {s[3], s[4]}, s[5], floe.state().trans});
        } else {
            auto& state = floe.state();
            state.pos = {s[0], s[1]};
            state.theta = s[2];
            state.speed = {s[3], s[4]};
            state.rot = s[5];
        }
        floe.reset_impulse(s[6]);
        m_states_origin[iter.first] = msg.mpi_source();
    }
}


}} // namespace floe::floes


#endif // FLOES_PARTIAL_FLOE_GROUP_HPP
