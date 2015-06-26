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
    typename TFloe,
    typename TFloeGroup_h = floe::variable::FloeGroup_h<>
>
class FloeGroup
{

public:

    //! Default constructor.
    // FloeGroup(){

    std::vector<TFloe> m_list_floe;

    inline void load_matlab_config(std::string filename) {
        using namespace floe::io::matlab;
        MatlabListSolid<double> list_so;
        cout << "Reading \"" << filename << "\" ... " << endl;
        read_list_so_from_file( filename, list_so);

        cout << "Importing floes ... " << endl;
        // m_list_floe = list_so_to_floes<TFloe>( list_so );
        for ( auto& floe : list_so_to_floes<TFloe>( list_so ) )
            add_floe(floe);
    };

    inline void add_floe( TFloe& floe )
        {
            m_list_floe.push_back(std::move(floe));
            m_floe_group_h.add_floe(&floe.get_floe_h());
        }

    inline TFloeGroup_h const& get_floe_group_h() const { return m_floe_group_h; }

private:

    TFloeGroup_h m_floe_group_h;

};

}} // namespace floe::variable


#endif // VARIABLE_FLOES_HPP
