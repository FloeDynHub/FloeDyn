/*!
 * \file floe/floes/identifiable_mixin.hpp
 * \brief Mixin adding identifier to any class
 * \author Roland Denis
 */

#ifndef FLOE_FLOES_IDENTIFIABLE_MIXIN_HPP
#define FLOE_FLOES_IDENTIFIABLE_MIXIN_HPP

namespace floe { namespace floes
{

/*! Mixin adding identifier to any class
 *
 * This adds id accessors for any class without hidding his own methods and properties.
 *
 * Example:
 * Identifiable<size_t, KinematicFloe> floe;
 * floe.id() = 3;
 * floe.attach_geometr_ptr(geo);
 * ...
 *
 * \tparam TIdentity Type of id (eg std::size_t)
 * \tparma TOwner    Type of identifiable class
 */
template <
    typename TIdentity,
    typename TOwner
>
class Identifiable 
    : public TOwner
{

public:
    //! TOwner constructors
    using TOwner::TOwner;

    //! Identity accessors
    TIdentity const&    get_id()                        const   { return m_id; }
    void                set_id( TIdentity const& id )           { m_id = id; }
    TIdentity const&    id()                            const   { return m_id; }
    TIdentity &         id()                                    { return m_id; }

private:
    TIdentity m_id;
};

}} // namespace floe::floes

#endif // FLOE_FLOES_IDENTIFIABLE_MIXIN_HPP

