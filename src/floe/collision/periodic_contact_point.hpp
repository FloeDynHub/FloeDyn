/*!
 * \file floe/collision/periodic_contact_point.hpp
 * \brief Contact point between two floes, may be via ghost floe in the case of periodic space
 * \author Quentin Jouet
 */

#ifndef FLOE_COLLISION_PERIODIC_CONTACT_POINT_HPP
#define FLOE_COLLISION_PERIODIC_CONTACT_POINT_HPP

#include "floe/collision/contact_point.hpp"

namespace floe { namespace collision
{

/*! Contact point between two floes, 
 * extended to the case of contact with ghost floes (periodic boundary conditions)
 *
 * \tparam TFloe        Floe type.
 * \tparam TGhostFloe   Ghost Floe type
 * \tparam TPoint       Point type.
 * \tparam TFrame       Frame type.
 */
 
template <
    typename TFloe,
    typename TGhostFloe,
    typename TPoint = typename TFloe::point_type,
    typename TFrame = floe::geometry::frame::UVFrame<TPoint>
>
struct PeriodicContactPoint : public ContactPoint<TFloe>
{

    using base_class = ContactPoint<TFloe>;

    //! Default constructor
    PeriodicContactPoint() : base_class() {}

    PeriodicContactPoint( TFloe const* floe1, TFloe const* floe2, TPoint const& pt1, TPoint const& pt2 ) : 
        base_class(floe1, floe2, pt1, pt2), floe1_trans{0, 0}, floe2_trans{0, 0} {}

    PeriodicContactPoint( TFloe const* floe1, TGhostFloe const* floe2, TPoint const& pt1, TPoint const& pt2 ) : 
        base_class(floe1, &floe2->original(), pt1, pt2), floe1_trans{0, 0}, floe2_trans{floe2->translation()} {}

    PeriodicContactPoint( TGhostFloe const* floe1, TFloe const* floe2, TPoint const& pt1, TPoint const& pt2 ) : 
        base_class(&floe1->original(), floe2, pt1, pt2), floe1_trans{floe1->translation()}, floe2_trans{0, 0} {}

    TPoint r1() const { 
        return base_class::frame.center() - base_class::floe1->state().pos - floe1_trans;
    }
    TPoint r2() const { 
        return base_class::frame.center() - base_class::floe2->state().pos - floe2_trans;
    }

    TPoint floe1_trans; //!< floe1 translation (translation from original of the ghost causing this contact)
    TPoint floe2_trans; //!< floe2 translation
};


}} // namespace floe::collision

#endif // FLOE_COLLISION_PERIODIC_CONTACT_POINT_HPP
