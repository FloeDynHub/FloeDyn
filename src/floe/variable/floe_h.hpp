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
#include "floe/variable/floe_alg.hpp"


namespace floe { namespace variable
{

/*! Floe_h
 *
 * It represents a single discret Floe.
 *
 */
template <
    typename TMesh     = geometry::TriangleMesh<geometry::Point<double>>,
    typename TFloe_alg = floe::variable::Floe_alg
>
class Floe_h
{

public:

    using floe_alg_type = TFloe_alg;

    // inline TMesh& kinematic_mesh() const { return m_kinematic_mesh; }
    TMesh m_kinematic_mesh;
    TMesh m_static_mesh;

    //! Floe_alg accessor
    inline floe_alg_type& get_floe_alg() { return m_floe_alg; }

private:

    floe_alg_type m_floe_alg;

};

}} // namespace floe::variable


#endif // VARIABLE_FLOE_H_HPP
