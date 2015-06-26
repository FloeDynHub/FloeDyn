/*!
 * \file variable/floe_group_h.hpp
 * \brief Discret Floe Configuration class
 * \author Quentin Jouet
 */

#ifndef VARIABLE_FLOES_H_HPP
#define VARIABLE_FLOES_H_HPP


#include "floe/variable/floe_group_alg.hpp"
#include "floe/variable/floe_h.hpp"


namespace floe { namespace variable
{

/*! FloeGroup_h
 *
 * It represents the floe configuration at the discrete level.
 *
 */
template <
    typename TFloe_h = floe::variable::Floe_h<>,
    typename TFloeGroup_alg = floe::variable::FloeGroup_alg<>
>
class FloeGroup_h
{

public:

    //! Default constructor.
    // FloeGroup_h()

    std::vector<TFloe_h*> m_list_floe_h;

    inline void add_floe( TFloe_h* floe )
        {
            m_list_floe_h.push_back(floe);
            m_floe_group_alg.add_floe(&floe->get_floe_alg());
        }

private:

    TFloeGroup_alg m_floe_group_alg;


};

}} // namespace floe::variable


#endif // VARIABLE_FLOES_H_HPP
