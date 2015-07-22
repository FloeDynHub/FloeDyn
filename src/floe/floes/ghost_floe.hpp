/*!
 * \file floe/floes/kinematic_floe.hpp
 * \brief 
 * \author Quentin Jouet
 */

#ifndef FLOE_FLOES_GHOST_FLOE_HPP
#define FLOE_FLOES_GHOST_FLOE_HPP


#include "floe/geometry/geometry.hpp"
#include "floe/geometry/arithmetic/point_operators.hpp"
#include "floe/geometry/frame/theta_frame.hpp"   
#include "floe/floes/floe_interface.hpp"

namespace floe { namespace floes
{

/*! Reflection of a floe (in periodic space)
 *
 */
template <
    typename TFloe
>
class GhostFloe : public FloeInterface<typename TFloe::static_floe_type, typename TFloe::state_type>
{

public:

    using floe_type = TFloe;
    using point_type = typename floe_type::point_type;
    using value_type = typename floe_type::value_type;
    using frame_type = typename floe_type::frame_type;
    using geometry_type = typename floe_type::geometry_type;
    using state_type = typename floe_type::state_type;
    using floe_interface_type = FloeInterface<typename TFloe::static_floe_type, typename TFloe::state_type>;
    using translate_strategy_type = boost::geometry::strategy::transform::translate_transformer<value_type, 2,2>;

    //! Default constructor
    GhostFloe() = delete;

    //! Constructor from floe, translation and floe id
    GhostFloe( TFloe const& floe, point_type translation, std::size_t floe_id = 0 ) :
        m_original_id{floe_id},
        m_floe{&floe},
        m_translation{translation},
        m_translator{translate_strategy_type(m_translation.x, m_translation.y)}
        {};

    //! Frame accessor
    frame_type const& frame() const;
    //! Geometry accessor
    geometry_type const& geometry() const;
    //! State accessor
    state_type& state() const;

    floe_type const& original() const { return *m_floe; }
    point_type translation() const { return m_translation; }

    const std::size_t m_original_id; //! Original floe id in group

private:
    const floe_type* m_floe;   //! Original floe
    point_type m_translation; //! Translation compared to original floe
    translate_strategy_type m_translator; //! Translation strategy for geometry transformation

    mutable frame_type m_frame;
    mutable geometry_type m_geometry;
    mutable state_type m_state;
};


template<typename TFloe>
typename GhostFloe<TFloe>::frame_type const&
GhostFloe<TFloe>::frame() const
{ 
    m_frame = m_floe->frame();
    m_frame.center() += m_translation;
    return m_frame;
}

template<typename TFloe>
typename GhostFloe<TFloe>::geometry_type const&
GhostFloe<TFloe>::geometry() const 
{
    boost::geometry::transform(m_floe->geometry(), m_geometry, m_translator);
    return m_geometry;
}

template<typename TFloe>
typename GhostFloe<TFloe>::state_type&
GhostFloe<TFloe>::state() const 
{
    m_state = m_floe->state();
    m_state.pos += m_translation;
    return m_state;
}


}} // namespace floe::floes
#endif // FLOE_FLOES_GHOST_FLOE_HPP

