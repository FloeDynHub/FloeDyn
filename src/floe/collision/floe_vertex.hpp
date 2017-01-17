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
	FloeVertex() : floe{nullptr} {}
	FloeVertex(TFloe const* floe_ptr) : floe{floe_ptr} {}
	TFloe const* floe;
	//! vertex id in parent graph (used with subgraphs)
    mutable std::size_t parent_descriptor;
};


}} // namespace floe::collision

#endif // FLOE_COLLISION_FLOE_VERTEX_HPP
