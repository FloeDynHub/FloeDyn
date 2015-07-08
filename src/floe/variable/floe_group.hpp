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
#include "floe/io/hdf5_writer.hpp"

#include "H5Cpp.h"
#include <algorithm>
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif


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
    using out_manager_type = io::HDF5Writer<FloeGroup<TFloe>>;

    //! Default constructor.
    // FloeGroup() : {}


    void load_matlab_config(std::string filename);

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
    value_type kinetic_energy();

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


template<typename TFloe>
typename FloeGroup<TFloe>::value_type
FloeGroup<TFloe>::kinetic_energy()
{
    return std::accumulate(
        m_list_floe.begin(), m_list_floe.end(), 0. , 
        [](value_type partial_sum, floe_type& floe) {return partial_sum + floe.kinetic_energy(); }
    );
}


}} // namespace floe::variable


#endif // VARIABLE_FLOES_HPP
