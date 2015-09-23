/*!
 * \file variable/floe_group.hpp
 * \brief Floe Configuration class
 * \author Quentin Jouet
 */

#ifndef VARIABLE_FLOES_HPP
#define VARIABLE_FLOES_HPP

#include <iostream>
#include "floe/variable/floe_group_h.hpp"

#include "floe/io/matlab/list_so_to_floes.hpp"
#include "floe/io/matlab/list_so_import.hpp"
#include "floe/io/matlab/list_so.hpp"
// #include "floe/io/false_hdf5_manager.hpp" // for gcc/MacOS
#include "floe/io/hdf5_manager.hpp"

namespace floe { namespace variable
{

/*! FloeGroup
 *
 * It represents the set of Floes.
 *
 */


template <
    typename TFloe
>
class FloeGroup
{

public:

    using floe_type = TFloe;
    using floe_group_h_type = floe::variable::FloeGroup_h<
        typename floe_type::floe_h_type
    >;
    using value_type = typename TFloe::value_type;
    using point_type = typename TFloe::point_type;
    using out_manager_type = io::HDF5Manager<FloeGroup<TFloe>>;

    //! Default constructor.
    // FloeGroup() : {}


    void load_matlab_config(std::string filename);

    double recover_states_from_file(std::string filename, double t);

    void out_hdf5(value_type time);

    inline void add_floe( floe_type& floe )
        {
            m_list_floe.push_back(std::move(floe));
            // auto& f = m_list_floe.back(); // TODO Why is this different from *1 ?
            // f.update();
            // m_floe_group_h.add_floe(f.get_floe_h());
        }

    inline floe_group_h_type const& get_floe_group_h() const { return m_floe_group_h; }
    inline floe_group_h_type& get_floe_group_h() { return m_floe_group_h; }
    inline std::vector<floe_type> const& get_floes() const { return m_list_floe; }
    inline std::vector<floe_type>& get_floes() { return m_list_floe; }

    //! kinetic energy of the group
    value_type kinetic_energy() const;
    //! sum all floes areas
    value_type total_area() const;
    //! sum all floes masses
    value_type total_mass() const;
    //! mass center of the group
    point_type mass_center() const;

private:

    std::vector<floe_type> m_list_floe;
    floe_group_h_type m_floe_group_h;
    out_manager_type m_out_manager;

};


template <
    typename TFloe
>
void FloeGroup<TFloe>::load_matlab_config(std::string filename) {
    using namespace floe::io::matlab;
    MatlabListSolid<double> list_so;
    cout << "Reading \"" << filename << "\" ... " << endl;
    read_list_so_from_file( filename, list_so);

    cout << "Importing floes ... " << endl;
    // m_list_floe = list_so_to_floes<floe_type>( list_so );
    for ( auto& floe : list_so_to_floes<floe_type>( list_so ) )
        add_floe(floe);
    for ( auto& floe : m_list_floe ){ // *1 ; todo : avoid that step
        floe.update();
        m_floe_group_h.add_floe(floe.get_floe_h());
    }
};


template <
    typename TFloe
>
void FloeGroup<TFloe>::out_hdf5(value_type time) {
    m_out_manager.save_step(time, *this);
};

template <
    typename TFloe
>
double FloeGroup<TFloe>::recover_states_from_file(std::string filename, double t) {
    return m_out_manager.recover_states(filename, t, *this);
};


template<typename TFloe>
typename FloeGroup<TFloe>::value_type
FloeGroup<TFloe>::kinetic_energy() const
{
    return std::accumulate(
        m_list_floe.begin(), m_list_floe.end(), 0. , 
        [](value_type partial_sum, floe_type const& floe) { return partial_sum + floe.kinetic_energy(); }
    );
}

template<typename TFloe>
typename FloeGroup<TFloe>::value_type
FloeGroup<TFloe>::total_area() const
{
    return std::accumulate(
        m_list_floe.begin(), m_list_floe.end(), 0. , 
        [](value_type partial_sum, floe_type const& floe) { return partial_sum + floe.area(); }
    );
}

template<typename TFloe>
typename FloeGroup<TFloe>::value_type
FloeGroup<TFloe>::total_mass() const
{
    return std::accumulate(
        m_list_floe.begin(), m_list_floe.end(), 0. , 
        [](value_type partial_sum, floe_type const& floe) { return partial_sum + floe.mass(); }
    );
}

template<typename TFloe>
typename FloeGroup<TFloe>::point_type
FloeGroup<TFloe>::mass_center() const
{
    return std::accumulate(
        m_list_floe.begin(), m_list_floe.end(), point_type{0,0} , 
        [](point_type partial_sum, floe_type const& floe) { return partial_sum + floe.mass() * floe.state().real_position(); }
    ) / total_mass();
}


}} // namespace floe::variable


#endif // VARIABLE_FLOES_HPP
