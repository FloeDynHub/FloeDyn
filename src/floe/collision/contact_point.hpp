/*!
 * \file floe/collision/contact_point.hpp
 * \brief Contact point between two floes and associated functions.
 * \author Roland Denis
 */

#ifndef FLOE_COLLISION_CONTACT_POINT_HPP
#define FLOE_COLLISION_CONTACT_POINT_HPP

#include <cmath>
#include <utility>
// #include <memory>

#include "floe/geometry/frame/uv_frame.hpp"
#include "floe/geometry/core/access.hpp"
#include "floe/geometry/core/coordinate_type.hpp"

namespace floe { namespace collision
{

/*! Contact point between two floes.
 *
 * The second vector v of the frame is normal to the contact, 
 * th first vector u is then tangential, so that (u,v) is orthonormal and direct oriented.
 *
 * \tparam TFloe    Floe type.
 * \tparam TPoint   Point type.
 * \tparam TFrame   Frame type.
 */
template <
    typename TFloe,
    typename TPoint = typename TFloe::point_type,
    typename TFrame = floe::geometry::frame::UVFrame<TPoint>
>
struct ContactPoint
{
    typedef TPoint      point_type;
    typedef TFrame      frame_type;
    typedef TFloe       floe_type;
    using real_type = typename floe_type::real_type;

    //! Default constructor
    ContactPoint() : floe1{nullptr}, floe2{nullptr}, frame{}, dist{std::numeric_limits<real_type>::max()} {}

    /*! Constructor given the contact frame
     * 
     * \param floe1    pointer to the first floe.
     * \param floe2    pointer to the second floe.
     * \param frame contact frame.
     * \param dist   optional distance between floes at this point
     */
    ContactPoint( TFloe const* floe1, TFloe const* floe2, frame_type const& frame, real_type distance = 0 ) 
        : floe1{floe1}, floe2{floe2}, frame{frame}, dist{distance}
    {}

    /*! Constructor given the two contact points.
     * 
     * \param floe1    pointer to the first floe.
     * \param floe2    pointer to the second floe.
     * \param pt1   contact point on the first object.
     * \param pt2   contact point on the second object.
     */
    ContactPoint( TFloe const* floe1, TFloe const* floe2, point_type const& pt1, point_type const& pt2 )
        : floe1{floe1}, floe2{floe2}
    {
        using namespace floe::geometry;

        const point_type u { get<1>(pt2) - get<1>(pt1), get<0>(pt1) - get<0>(pt2) };
        typename coordinate_type<point_type>::type const norm_u = std::sqrt( std::pow(get<0>(u), 2) + std::pow(get<1>(u), 2) );

        frame = frame_type{ pt1, { get<0>(u)/norm_u, get<1>(u)/norm_u } };
        dist = norm_u;
    }

    //! Conversion to frame_type returning the frame
    inline explicit
    operator frame_type() const { return frame; }

    real_type relative_speed() const
    {
        if (*floe_states_changed) { calc_relative_speed(); *floe_states_changed = false; }
        return *m_relative_speed;
    }

    //! Calculate the relative speed of the two floes relatively to the contact frame
    void
    calc_relative_speed() const
    {
        using namespace floe::geometry;

        const auto& state1 = floe1->state();
        const auto& state2 = floe2->state();

        const point_type& speed1 = state1.speed + direct_orthogonal( r1() ) * state1.rot;
        const point_type& speed2 = state2.speed + direct_orthogonal( r2() ) * state2.rot;

        *m_relative_speed = dot_product(speed2 - speed1, frame.v());
    }

    virtual point_type r1() const { 
        return frame.center() - floe1->state().pos;
    }

    virtual point_type r2() const { 
        return frame.center() - floe2->state().pos;
    }


    //! Return true if the contact is active (relative speed of the two contact points is negative)
    inline bool is_active() const
    {
        // return relative_speed() < - dist / 50; // 10 = DT_DEFAULT // get dt_defaut ?
        return relative_speed() < 0;
    }

    TFloe const* floe1; //!< First floe in contact
    TFloe const* floe2; //!< Second floe in contact
    TFrame       frame; //!< Frame of contact
    real_type dist; //!< Distance between floes at this point
    mutable std::shared_ptr<real_type> m_relative_speed{std::make_shared<real_type>(-1)}; //!< Relative speed cash (for performances)
    //! Floe states changed (Need new relative speed calculation), shared pointer to be shared with subgraphs.
    mutable std::shared_ptr<bool> floe_states_changed{std::make_shared<bool> (true)};
};

//! Exterior function to test if a contact is active (should be the only place ...)
template < typename TContactPoint >
inline bool is_active( TContactPoint const& contact )
{
    return contact.is_active();
}

}} // namespace floe::collision

#endif // FLOE_COLLISION_CONTACT_POINT_HPP
