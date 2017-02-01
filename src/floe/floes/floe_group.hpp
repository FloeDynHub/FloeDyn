/*!
 * \file floes/floe_group.hpp
 * \brief Floe Configuration class
 * \author Quentin Jouet
 */

#ifndef VARIABLE_FLOES_HPP
#define VARIABLE_FLOES_HPP

#include <iostream>
#include "floe/floes/floe_group_h.hpp"
#include "floe/io/hdf5_config_import.hpp"

#include "floe/io/matlab/list_so_to_floes.hpp"
#include "floe/io/matlab/list_so_import.hpp"
#include "floe/io/matlab/list_so.hpp"

// #include "floe/io/inter_process_message.hpp" // TODO move out of non-MPI code ?

#include <algorithm>
#include <random>

namespace floe { namespace floes
{

/*! FloeGroup
 *
 * It represents the set of Floes.
 *
 */
template <
    typename TFloe,
    typename TFloeList = std::vector<TFloe>
>
class FloeGroup
{

public:

    using floe_type = TFloe;
    using floe_group_h_type = floe::floes::FloeGroup_h<
        typename floe_type::floe_h_type
    >;
    using value_type = typename TFloe::value_type;
    using point_type = typename TFloe::point_type;
    // using floe_list_type = std::vector<floe_type>;
    using floe_list_type = TFloeList;
    using window_type = std::array<value_type, 4>;

    //! Default constructor.
    FloeGroup() : m_window{{0,0,0,0}} {}

    //! Load floes and initial states from matlab file
    void load_matlab_config(std::string filename);
    //! Load floes and initial states from hdf5 file
    void load_h5_config(std::string filename);
    virtual void post_load_floe(){;}

    // Accessors
    inline floe_group_h_type const& get_floe_group_h() const { return m_floe_group_h; }
    inline floe_group_h_type& get_floe_group_h() { return m_floe_group_h; }
    // inline floe_list_type const& get_all_floes() const { return m_list_floe; }
    // inline floe_list_type& get_all_floes() { return m_list_floe; }
    // virtual inline floe_list_type const& get_floes() const { return get_all_floes(); }
    // virtual inline floe_list_type& get_floes() { return get_all_floes(); }
    virtual inline floe_list_type const& get_floes() const { return m_list_floe; }
    virtual inline floe_list_type& get_floes() { return m_list_floe; }
    inline std::size_t nb_floes() { return this->get_floes().size(); }
    inline window_type const& get_initial_window() const { return m_window; }
    void set_initial_window(window_type win){ m_window = win; }
    value_type ocean_window_area() const;
    value_type initial_window_area() const { return (m_window[1] - m_window[0]) * (m_window[3] - m_window[2]); }

    //! kinetic energy of the group
    value_type kinetic_energy() const;
    //! sum all floe areas
    value_type total_area() const;
    //! sum all floe masses
    value_type total_mass() const;
    //! mass center of the group
    point_type mass_center() const;
    //! bounding window of floe group (return array of min_x, max_x, min_y, max_y)
    std::array<value_type, 4> bounding_window(value_type margin = 1) const;
    //! bounding window area
    value_type bounding_window_area(value_type margin = 1) const;
    //! Floe Concentration
    value_type floe_concentration() const;
    //! Initial Floe Concentration
    value_type initial_concentration() const;

    //! Stops floes contained in centered window (for generator)
    void stop_floes_in_window(value_type width, value_type height);
    //! Change all floes friction coeff (generator VS simulation)
    void set_mu_static(value_type mu_static);
    //! Reset received impulse
    void reset_impulses() const;

    //! save floe states for prospective recover
    void backup_step_states();
    //! Load saved previous time step floe states
    virtual void recover_previous_step_states();
    //! Set floes thickness to normal distributed random values around default
    void randomize_floes_thickness(value_type coeff);
    //! Set floes oceanic skin drag coeff to normal distributed random values around default
    void randomize_floes_oceanic_skin_drag(value_type coeff);

    // using message_type = io::InterProcessMessage<value_type>;
    //! For MPI, states received from another process via inter-process-message
    // virtual void update_floe_states(message_type const& msg, bool update=true);

    virtual inline int absolute_id(int id) const { return id; }
    
protected:

    floe_list_type m_list_floe; //!< List of floes
    floe_group_h_type m_floe_group_h; //!< Discrete floe group (access to floes discretisation)
    //! initial reference window (min_x, max_x, min_y, max_y)
    window_type m_window;
    //! floes previous states backup
    std::vector<typename floe_type::state_type> m_previous_step_states;

};


template <typename TFloe, typename TFloeList>
void FloeGroup<TFloe, TFloeList>::load_matlab_config(std::string filename) {
    using namespace floe::io::matlab;
    MatlabListSolid<double> list_so;
    cout << "Reading \"" << filename << "\" ... " << endl;
    read_list_so_from_file( filename, list_so);
    cout << "Importing floes ... " << endl;
    list_so_to_floes( list_so, m_list_floe );
    for ( auto& floe : m_list_floe ){
        floe.update();
        m_floe_group_h.add_floe(floe.get_floe_h());
    }
    this->post_load_floe();
};

template <typename TFloe, typename TFloeList>
void FloeGroup<TFloe, TFloeList>::load_h5_config(std::string filename) {
    std::cout << "Reading \"" << filename << "\" ... " << std::endl;
    floe::io::import_floes_from_hdf5(filename, *this);
    this->post_load_floe();
};

template <typename TFloe, typename TFloeList>
typename FloeGroup<TFloe, TFloeList>::value_type
FloeGroup<TFloe, TFloeList>::kinetic_energy() const
{
    return std::accumulate(
        get_floes().begin(), get_floes().end(), 0. , 
        [](value_type partial_sum, floe_type const& floe) { return partial_sum + floe.kinetic_energy(); }
    );
}

template <typename TFloe, typename TFloeList>
typename FloeGroup<TFloe, TFloeList>::value_type
FloeGroup<TFloe, TFloeList>::total_area() const
{
    return std::accumulate(
        get_floes().begin(), get_floes().end(), 0. , 
        [](value_type partial_sum, floe_type const& floe) { return partial_sum + floe.area(); }
    );
}

template <typename TFloe, typename TFloeList>
typename FloeGroup<TFloe, TFloeList>::value_type
FloeGroup<TFloe, TFloeList>::total_mass() const
{
    return std::accumulate(
        get_floes().begin(), get_floes().end(), 0. , 
        [](value_type partial_sum, floe_type const& floe) { return partial_sum + floe.mass(); }
    );
}

template <typename TFloe, typename TFloeList>
typename FloeGroup<TFloe, TFloeList>::point_type
FloeGroup<TFloe, TFloeList>::mass_center() const
{
    return std::accumulate(
        get_floes().begin(), get_floes().end(), point_type{0,0} , 
        [](point_type partial_sum, floe_type const& floe) { return partial_sum + floe.mass() * floe.state().real_position(); }
    ) / total_mass();
}

template <typename TFloe, typename TFloeList>
std::array<typename FloeGroup<TFloe, TFloeList>::value_type, 4>
FloeGroup<TFloe, TFloeList>::bounding_window(value_type margin) const
{
    value_type min_x, min_y, max_x, max_y;
    min_x = min_y = std::numeric_limits<value_type>::max();
    max_x = max_y = - std::numeric_limits<value_type>::max();

    for (auto const& floe : get_floes())  
        for (auto const& pt : floe.geometry().outer())
        {
            min_x = std::min(min_x, pt.x);
            min_y = std::min(min_y, pt.y);
            max_x = std::max(max_x, pt.x);
            max_y = std::max(max_y, pt.y);
        }
    return {{min_x - margin, max_x + margin, min_y - margin, max_y + margin}};
}

template <typename TFloe, typename TFloeList>
typename FloeGroup<TFloe, TFloeList>::value_type
FloeGroup<TFloe, TFloeList>::bounding_window_area(value_type margin) const {
    auto a = bounding_window(margin);
    return ((a[1] - a[0]) * (a[3] - a[2]));
}

template <typename TFloe, typename TFloeList>
typename FloeGroup<TFloe, TFloeList>::value_type
FloeGroup<TFloe, TFloeList>::floe_concentration() const {
    return total_area() / bounding_window_area(0);
}

template <typename TFloe, typename TFloeList>
typename FloeGroup<TFloe, TFloeList>::value_type
FloeGroup<TFloe, TFloeList>::initial_concentration() const {
    return total_area() / initial_window_area();
}

template <typename TFloe, typename TFloeList>
void FloeGroup<TFloe, TFloeList>::stop_floes_in_window(value_type width, value_type height)
{
    for (auto& floe : get_floes())
    {
        auto const& out = floe.geometry().outer();
        if (std::all_of(out.begin(), out.end(), [&](point_type const& pt){ return (std::abs(pt.x) < width / 2 && std::abs(pt.y) < height / 2); }))
        {
            floe.state().speed = point_type{0,0}; floe.state().rot = 0;
        }
    }   
}

template <typename TFloe, typename TFloeList>
void FloeGroup<TFloe, TFloeList>::set_mu_static(value_type mu_static)
{
    for (auto& floe : get_floes()) floe.set_mu_static(mu_static);
}

template <typename TFloe, typename TFloeList>
void FloeGroup<TFloe, TFloeList>::reset_impulses() const
{
    for (auto& floe : get_floes()) floe.reset_impulse();
}

template <typename TFloe, typename TFloeList>
typename FloeGroup<TFloe, TFloeList>::value_type
FloeGroup<TFloe, TFloeList>::ocean_window_area() const
{
    if (initial_window_area()) return initial_window_area();
    else return bounding_window_area();
}

template <typename TFloe, typename TFloeList>
void
FloeGroup<TFloe, TFloeList>::backup_step_states()
{
    m_previous_step_states.resize(m_list_floe.size());
    for (std::size_t i = 0; i < m_list_floe.size(); ++i)
    {
        m_previous_step_states[i] = m_list_floe[i].state();
    }
}


template <typename TFloe, typename TFloeList>
void
FloeGroup<TFloe, TFloeList>::recover_previous_step_states()
{
    for (std::size_t i = 0; i < m_list_floe.size(); ++i)
    {
        m_list_floe[i].set_state(m_previous_step_states[i]);
    }
}

template <typename TFloe, typename TFloeList>
void
FloeGroup<TFloe, TFloeList>::randomize_floes_thickness(value_type coeff)
{
    auto dist = std::normal_distribution<value_type>{1, coeff};
    auto gen = std::default_random_engine{};
    for (auto& floe : get_floes()){
        auto& static_floe = floe.static_floe();
        static_floe.set_thickness(static_floe.thickness() * dist(gen));
    }
}

template <typename TFloe, typename TFloeList>
void
FloeGroup<TFloe, TFloeList>::randomize_floes_oceanic_skin_drag(value_type coeff)
{
    auto dist = std::normal_distribution<value_type>{1, coeff};
    auto gen = std::default_random_engine{};
    gen.seed(1); // to get a different set of values than for floes thickness
    for (auto& floe : get_floes()){
        auto& static_floe = floe.static_floe();
        static_floe.set_C_w(static_floe.C_w() * dist(gen));
    }
}

}} // namespace floe::floes


#endif // VARIABLE_FLOES_HPP