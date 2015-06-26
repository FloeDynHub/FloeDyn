/*!
 * \file variable/floe.hpp
 * \brief Floe
 * \author Quentin Jouet

 USELESS FILE (use static/kinematic floe instead)

 */

#ifndef VARIABLE_FLOE_HPP
#define VARIABLE_FLOE_HPP


#include "floe/floes/static_floe.hpp"


namespace floe { namespace variable
{

/*! Floe
 *
 * It represents a single Floe.
 * USELESS
 *
 */
template <
    typename TStaticFloe = floe::floes::StaticFloe,
    typename TKinematicFloe = floe::floes::KinematicFloe<TStaticFloe>,
>
class Floe
{

public:

    //! Default constructor.
    // Floes()

    TKinematicFloe m_floe;

private:

    TFloe_h m_floes_h;

};

}} // namespace floe::variable


#endif // VARIABLE_FLOE_HPP
