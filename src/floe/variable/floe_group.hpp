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

    //! Default constructor.
    // FloeGroup() : {}

    void load_matlab_config(std::string filename);

    inline floe_group_h_type const& get_floe_group_h() const { return m_floe_group_h; }
    inline floe_group_h_type& get_floe_group_h() { return m_floe_group_h; }
    inline std::vector<floe_type> const& get_floes() const { return m_list_floe; }
    inline std::vector<floe_type>& get_floes() { return m_list_floe; }

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
    

private:

    std::vector<floe_type> m_list_floe;
    floe_group_h_type m_floe_group_h;

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
    list_so_to_floes( list_so, m_list_floe );
    for ( auto& floe : m_list_floe ){
        floe.update();
        m_floe_group_h.add_floe(floe.get_floe_h());
    }
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

template<typename TFloe>
std::array<typename FloeGroup<TFloe>::value_type, 4>
FloeGroup<TFloe>::bounding_window(value_type margin) const
{
    value_type min_x, min_y, max_x, max_y;
    min_x = min_y = std::numeric_limits<value_type>::max();
    max_x = max_y = - std::numeric_limits<value_type>::max();

    for (auto const& floe : m_list_floe)  
        for (auto const& pt : floe.geometry().outer())
        {
            min_x = std::min(min_x, pt.x);
            min_y = std::min(min_y, pt.y);
            max_x = std::max(max_x, pt.x);
            max_y = std::max(max_y, pt.y);
        }
    return {{min_x - margin, max_x + margin, min_y - margin, max_y + margin}};
}


}} // namespace floe::variable


#endif // VARIABLE_FLOES_HPP
