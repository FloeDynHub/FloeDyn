/*!
 * \file floe/floes/floe_interface.hpp
 * \brief Virtual interface of a floe.
 * \author Quentin Jouet
 */

#ifndef FLOE_FLOES_FLOE_INTERFACE_HPP
#define FLOE_FLOES_FLOE_INTERFACE_HPP

namespace floe { namespace floes
{

/*! public interface of a floe
 */
template< typename TFloe, typename TState = typename TFloe::state_type >
class FloeInterface
{

public:
    using floe_type = TFloe;
    using frame_type = typename floe_type::frame_type;
    using geometry_type = typename floe_type::geometry_type;
    using state_type = TState;

    //! Frame accessor
    virtual frame_type const&  frame() const = 0;
    //! Geometry accessor
    virtual geometry_type const& geometry() const = 0;
    //! State accessor
    virtual state_type const& state() const = 0;

};

}} // namespace floe::floes
#endif // FLOE_FLOES_FLOE_INTERFACE_HPP

