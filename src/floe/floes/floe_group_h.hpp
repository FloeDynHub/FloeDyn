/*!
 * \file floes/floe_group_h.hpp
 * \brief Discret Floe Configuration class
 * \author Quentin Jouet
 */

#ifndef VARIABLE_FLOES_H_HPP
#define VARIABLE_FLOES_H_HPP


#include "floe/floes/floe_h.hpp"


namespace floe { namespace floes
{

/*! FloeGroup_h
 *
 * It represents the floe configuration at the discrete level.
 *
 */
template <
    typename TFloe_h
>
class FloeGroup_h
{

public:

    //! Default constructor.
    // FloeGroup_h()

    std::vector<TFloe_h*> m_list_floe_h; //!< Discrete floes list

    //! Add a floe in the list
    inline void add_floe( TFloe_h& floe )
        {
            m_list_floe_h.push_back(&floe);
        }

};

}} // namespace floe::floes


#endif // VARIABLE_FLOES_H_HPP
