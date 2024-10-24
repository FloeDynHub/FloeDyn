/*!
 * \file floe/floes/partial_floe_group.hpp
 * \brief Partial Floe Configuration class (subset of floes)
 * \author Quentin Jouet
 */

#ifndef FLOES_PARTIAL_FLOE_GROUP_HPP
#define FLOES_PARTIAL_FLOE_GROUP_HPP

#define WHEREAMI std::cout << std::endl << "no crash until line " << __LINE__ << " in the file " __FILE__ << std::endl;



#include "floe/floes/floe_group.hpp"
#include "floe/arithmetic/filtered_container.hpp"
#include "floe/io/inter_process_message.hpp"
#include "floe/generator/mesh_generator.hpp"

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
    using geometry_type = typename floe_type::geometry_type;
    using point_type = typename floe_type::point_type;
    using static_floe_type = typename floe_type::static_floe_type;
    using mesh_type = typename floe_type::mesh_type;
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

    // fracture !
    void add_floe(geometry_type geometry, std::size_t parent_floe_idx);
    void fracture_biggest_floe();
    // size_t fracture_above_threshold(real_type threshold);
    int fracture_floes();
    void melt_floes();
    void update_list_ids_active();//{std::cout<<"test"<<std::endl;}

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

template <typename TFloe, typename TFloeList>
void
PartialFloeGroup<TFloe, TFloeList>::fracture_biggest_floe()
{
	real_type max_area = 0;
    int biggest_floe_idx = 0;
	for (std::size_t i = 0; i < base_class::get_floes().size(); ++i){
        auto& floe = base_class::get_floes()[i];
        if (!floe.is_obstacle() && floe.area() > max_area){
            max_area = floe.area();
            biggest_floe_idx = i;
        }
    }

    auto new_geometries = base_class::get_floes()[biggest_floe_idx].fracture_floe();
    for (std::size_t i = 0; i < new_geometries.size(); ++i){
    	this->add_floe(new_geometries[i], biggest_floe_idx);
    }

    // Desactivate cracked floe
    base_class::get_floes()[biggest_floe_idx].state().desactivate();

    this->update_list_ids_active();

    for (auto & floe : this->get_floes()) { // TODO why is it needed ?
        floe.static_floe().attach_mesh_ptr(&floe.get_floe_h().m_static_mesh);
        floe.update();
    }
}

// template <typename TFloe, typename TFloeList>
// size_t
// PartialFloeGroup<TFloe, TFloeList>::fracture_above_threshold(real_type threshold)
// {
// 	// returns the number of cracked floes ``
//     size_t nCracked(0);
// 	for (std::size_t iFloe = 0; iFloe < base_class::get_floes().size(); ++iFloe){
//         auto& floe = base_class::get_floes()[iFloe];
//         // WHEREAMI
//         if (!floe.prepare_elasticity())
//             std::cout << "FEM computation initialization failed" << std::endl;
//         // WHEREAMI
//         if (!floe.is_obstacle() && floe.total_received_impulse() > 0)
//         {
//             // WHEREAMI
//             std::cout << "trying to solve elasticity " << std::endl;
//             if (!floe.solve_elasticity())
//                 std::cout << "Solve on floe " << iFloe << " has failed." << std::endl;
//             // WHEREAMI
//         }
//         if (!floe.is_obstacle() && floe.total_received_impulse() > threshold){
//             // if the impulse is greater than a threshold, flow is fractured
//             auto new_geometries = floe.fracture_floe();
//             // WHEREAMI
//             std::cout << "fracturing floe " << iFloe << " whose impulse reaches " << floe.total_received_impulse()<< " replaced by " << new_geometries.size() << " new geometries" << std::endl;
//             for (std::size_t i = 0; i < new_geometries.size(); ++i){
//                 // new geometries are added to the floe list.
//                 // note : iFloe is needed to initialize correctly the new floes states
//                 this->add_floe(new_geometries[i], iFloe);
//                 // WHEREAMI
//             }
//             // WHEREAMI
//             // previous flow is deactivated
//             base_class::get_floes()[iFloe].state().desactivate(); // floe.state().desactivate(); does not work
//             nCracked++;
//             this->update_list_ids_active();

//             for (auto & floe : this->get_floes()) { // TODO why is it needed ?
//                 floe.static_floe().attach_mesh_ptr(&floe.get_floe_h().m_static_mesh);
//                 floe.update();
//             }
//             // WHEREAMI
//         }
//     }
//     return nCracked;
// }


template <typename TFloe, typename TFloeList>
int
PartialFloeGroup<TFloe, TFloeList>::fracture_floes()
{
    int n_fractured = 0;
    // real_type min_area(400);
    real_type min_area(0.001);
    std::map<std::size_t, std::vector<geometry_type>> all_new_geometries;
    for (std::size_t i = 0; i < base_class::get_floes().size(); ++i){
        auto& floe = base_class::get_floes()[i];
        if (floe.is_obstacle())
        {
            std::cout << "Ignoring Floe " << i << " (obstacle)." << std::endl;
            continue;
        }
        if (floe.area() < min_area)
        {
            std::cout << "Ignoring Floe " << i << " (too small). " << std::endl;
            continue;
        }
        if (!floe.has_been_impacted())
        {
            std::cout << "Ignoring Floe " << i << " (no impact). " << std::endl;
            continue;
        }
        auto new_geometries = base_class::get_floes()[i].fracture_floe_from_collisions();
        // auto new_geometries = floe.fracture_floe_from_collisions_fem();
        std::cout << "Looking for fracture in Floe " << i << ":";
        if (new_geometries.size() > 0){
            std::cout << " fractured in " << new_geometries.size() << " parts" << std::endl;
            all_new_geometries[i] = new_geometries;
            n_fractured++;
        }
        else{
            std::cout << " not fractured " << std::endl;
        }
    }

    // Add new floes
    for (auto const& iter : all_new_geometries){
        for (std::size_t j = 0; j < iter.second.size(); ++j){
            this->add_floe(iter.second[j], iter.first);
        }
    }
    // Desactivate cracked floe
    for (auto const& iter : all_new_geometries){
        base_class::get_floes()[iter.first].state().desactivate();
    }
    this->update_list_ids_active();

    // // Deactivate too small floes
    // for (auto & floe : base_class::get_floes()){
    //     if (floe.area() < min_area)
    //     {
    //         floe.state().desactivate();
    //         // std::cout << "Floe " << i << " is too small and has been deactivated." << std::endl;
    //         std::cout << "Floe is too small and has been deactivated." << std::endl;
    //     }
    // }
    // this->update_list_ids_active();

    for (auto & floe : this->get_floes()) { // TODO why is it needed ?
        floe.static_floe().attach_mesh_ptr(&floe.get_floe_h().m_static_mesh);
        floe.update();
        floe.unset_fem_problem_prepared();
    }
    return n_fractured;
}

template <typename TFloe, typename TFloeList>
void
PartialFloeGroup<TFloe, TFloeList>::melt_floes()
{
    // Dumb melting model for feature testing :
    // each floe looses random thickness ~ 0.3mm / time step
    auto dist = std::normal_distribution<real_type>{1, 0.1};
    auto gen = std::default_random_engine{};
	for (auto& floe : base_class::get_floes()){
        auto& static_floe = floe.static_floe();
        static_floe.set_thickness(static_floe.thickness() - 0.0003 * dist(gen));
        // Desactivate molten floes
        if (static_floe.thickness() < base_class::m_min_thickness) {
            std::cout << "MOLTEN !" << std::endl;
            static_floe.set_thickness(0);
            floe.state().desactivate();
        }
    }
    this->update_list_ids_active();
}

template <typename TFloe, typename TFloeList>
void
PartialFloeGroup<TFloe, TFloeList>::add_floe(geometry_type shape, std::size_t parent_floe_idx)
{
	// Resize floe group, set all floe properties
    auto& list_floes = base_class::get_floes();
    std::size_t parent_floe_abs_id = list_floes.absolute_id(parent_floe_idx);
    base_class::get_floes().filter_off();

    // Create Kinematic floe
    list_floes.resize(list_floes.size() + 1);
    auto& floe = list_floes[list_floes.size() - 1];

    // link static floe
    floe.attach_static_floe_ptr(std::unique_ptr<static_floe_type>(new static_floe_type()));
    auto& static_floe = floe.static_floe();

    // Create mesh
    auto mesh = floe::generator::generate_mesh_for_shape<geometry_type, mesh_type>(shape);

    // Center mesh and shape on new floe's center of mass
     using integration_strategy = floe::integration::RefGaussLegendre<real_type,2,2>;
     auto mass_center = floe::integration::integrate(
         [] (real_type x, real_type y) { return point_type{x, y}; },
         mesh, integration_strategy()
     ) / floe::integration::integrate(
         [] (real_type x, real_type y) { return 1.; },
         mesh,integration_strategy()
     );
    geometry_type shape_cpy = shape;
    geometry::transform( shape_cpy, shape, geometry::frame::transformer( typename floe_type::frame_type{-mass_center, 0} ));
    mesh_type mesh_cpy = mesh;
    geometry::transform( mesh_cpy, mesh, geometry::frame::transformer( typename floe_type::frame_type{-mass_center, 0} ));

    // Save mesh and shape
    std::unique_ptr<typename floe_type::geometry_type> geometry(new typename floe_type::geometry_type(shape));
    static_floe.attach_geometry_ptr(std::move(geometry));

    mesh_type& floe_mesh = floe.get_floe_h().m_static_mesh;
    floe_mesh = mesh;
    static_floe.attach_mesh_ptr(&floe_mesh);

    this->get_floe_group_h().add_floe(floe.get_floe_h());
    // Compute and set space-time state
    auto& parent_floe = list_floes[parent_floe_abs_id];
    point_type rotated_mc {
        mass_center.x * std::cos(parent_floe.state().theta) - mass_center.y * std::sin(parent_floe.state().theta),
        mass_center.x * std::sin(parent_floe.state().theta) + mass_center.y * std::cos(parent_floe.state().theta),
    };
    floe.set_state({
        parent_floe.state().pos + rotated_mc,
        parent_floe.state().theta,
        parent_floe.ice_speed(parent_floe.state().pos + rotated_mc),
        parent_floe.state().rot,
        parent_floe.state().trans
    });
    // FLoes's inherited caracteristics
    // Random floe oceanic skin drag variation
    auto dist = std::normal_distribution<real_type>{1, 0.02};
    auto gen = std::default_random_engine{};
    gen.seed(list_floes.absolute_size());
    static_floe.set_C_w(static_floe.C_w() * dist(gen));
    static_floe.set_thickness(parent_floe.static_floe().thickness() * dist(gen));
    base_class::get_floes().filter_on();
}

/*

template <typename TFloe, typename TFloeList>
real_type
PartialFloeGroup<TFloe, TFloeList>::max_floe_area()
{
	real_type max_area {0.0};
	for (std::size_t i = 0; i < base_class::get_floes().size(); ++i){
    	max_floe_are = std::max(max_area, base_class::get_floes()[i].static_floe.area());
    }
    return max_area;
}

*/

// method which checks which floe is active or not and creates a filter
// find a way to do that more efficiency
template <typename TFloe, typename TFloeList>
void
PartialFloeGroup<TFloe, TFloeList>::update_list_ids_active()
{
	base_class::get_floes().filter_off();
	std::vector<std::size_t> m_list_id_active_floe;
	// add active floe to the liste of indice of active floe
	for (std::size_t i = 0; i < base_class::get_floes().size(); ++i){
    	if ( base_class::get_floes()[i].state().is_active()) { m_list_id_active_floe.push_back(i); }
    }
    this->update_partial_list(m_list_id_active_floe);
    base_class::get_floes().filter_on();
}


}} // namespace floe::floes


#endif // FLOES_PARTIAL_FLOE_GROUP_HPP
