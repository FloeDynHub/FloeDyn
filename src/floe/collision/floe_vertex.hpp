/*!
 * \file floe/collision/floe_vertex.hpp
 * \brief Representation of a contact between two floes.
 * \author Quentin Jouet
 */

#ifndef FLOE_COLLISION_FLOE_VERTEX_HPP
#define FLOE_COLLISION_FLOE_VERTEX_HPP


namespace floe { namespace collision
{

/*! Contact between two floes.
 *
 * \tparam TContactPoint   contact point type.
 */

template <
    typename TFloe
>
struct FloeVertex
{
	using real_type = typename TFloe::real_type;
	FloeVertex() : floe{nullptr} {}
	FloeVertex(TFloe const* floe_ptr) : floe{floe_ptr} {}
	//! Get received impulse
	inline real_type impulse() const { return *m_impulse_received; }
	//! Set received impulse
	inline void add_impulse_received(real_type impulse) const { *m_impulse_received += impulse; }
	//! The floe
	TFloe const* floe;
	//! vertex id in parent graph (used with subgraphs)
	mutable std::size_t parent_descriptor;
	// mutable real_type m_impulse_received{0};
	mutable std::shared_ptr<real_type> m_impulse_received{std::make_shared<real_type>(0)};
};


}} // namespace floe::collision

#endif // FLOE_COLLISION_FLOE_VERTEX_HPP
