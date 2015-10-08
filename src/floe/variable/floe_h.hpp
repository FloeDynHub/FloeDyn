/*!
 * \file variable/floe_h.hpp
 * \brief Discret Floe
 * \author Quentin Jouet
 */

#ifndef VARIABLE_FLOE_H_HPP
#define VARIABLE_FLOE_H_HPP

#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/point.hpp"           // Default point type
#include "floe/geometry/geometries/triangle_mesh.hpp"   // Default mesh type


namespace floe { namespace variable
{

/*! Floe_h
 *
 * It represents a single discret Floe.
 *
 */
template <
    typename TMesh     = geometry::TriangleMesh<geometry::Point<double>>
>
class Floe_h
{

public:

    // inline TMesh& kinematic_mesh() const { return m_kinematic_mesh; }
    TMesh m_kinematic_mesh;
    TMesh m_static_mesh;


};

}} // namespace floe::variable


#endif // VARIABLE_FLOE_H_HPP
