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
    typename TFloe_h
>
class FloeGroup_h
{

public:

    using floe_group_alg_type = floe::variable::FloeGroup_alg<
        typename TFloe_h::floe_alg_type
    >;

    //! Default constructor.
    // FloeGroup_h()

    std::vector<TFloe_h*> m_list_floe_h;

    inline void add_floe( TFloe_h& floe )
        {
            m_list_floe_h.push_back(&floe);
            m_floe_group_alg.add_floe(floe.get_floe_alg());
        }

    inline floe_group_alg_type const& get_floe_group_alg() const { return m_floe_group_alg; }
    inline floe_group_alg_type& get_floe_group_alg() { return m_floe_group_alg; }

private:

    floe_group_alg_type m_floe_group_alg;



};

}} // namespace floe::variable


#endif // VARIABLE_FLOES_H_HPP
