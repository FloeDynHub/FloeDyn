/*!
 * \file variable/floe_group_alg.hpp
 * \brief Algebraic Floe Configuration class
 * \author Quentin Jouet
 */

#ifndef VARIABLE_FLOES_ALG_HPP
#define VARIABLE_FLOES_ALG_HPP


#include "floe/variable/floe_alg.hpp"


namespace floe { namespace variable
{

/*! FloeGroup_alg
 *
 * It represents the floe configuration at the algebraic level.
 *
 */
template <
    typename TFloe_alg = floe::variable::Floe_alg
>
class FloeGroup_alg
{

public:

    //! Default constructor.
    // FloeGroup_alg()

    std::vector<TFloe_alg*> m_list_floe_alg;

    inline void add_floe( TFloe_alg* floe )
        {
        }

private:

};

}} // namespace floe::variable


#endif // VARIABLE_FLOES_ALG_HPP
