/*!
 * \file floe/io/matlab/list_so_to_floes.hpp
 * \brief Translation from the list_so structure to a list of kinematic floes.
 * \author Roland Denis
 */

#ifndef FLOE_IO_MATLAB_LIST_SO_TO_FLOES_HPP
#define FLOE_IO_MATLAB_LIST_SO_TO_FLOES_HPP

#include <vector>
#include <stdexcept>

#include "floe/io/matlab/list_so.hpp"
#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"

namespace floe { namespace io { namespace matlab
{

/*! Convert a MatlabListSolid structure (e.g. list_so) to a list of Floes.
 *
 * \tparam TKinematicFloe Type of the kinematic floe.
 * \tparam TMatlabList    Type of the Matlab list_so structure.
 * \param  list_so        The list_so structure.
 *
 * \todo Instead having a template for list_so, use only a template for T and accept MatlabListSolid<T> parameter.
 */
template < 
    typename TKinematicFloe,
    typename TMatlabListSolid
>
std::vector<TKinematicFloe*>
list_so_to_floes( TMatlabListSolid const& list_so )
{
    std::vector<TKinematicFloe*> floes;

    // Typedefs
    using TStaticFloe = typename TKinematicFloe::static_floe_type;
    using real        = typename TKinematicFloe::value_type;
    using TGeometry   = typename TKinematicFloe::geometry_type;
    using TMesh       = typename TKinematicFloe::mesh_type;
    using TState      = typename TKinematicFloe::state_type;

    // Verify that there is the same number of Solid and Movement object
    if ( list_so.geo.size() != list_so.mov.size() )
        throw std::length_error("geo and mov list are of different sizes.");

    std::size_t n_floes = list_so.geo.size();
    floes.reserve(n_floes);

    // Import each floe
    for ( std::size_t i = 0; i < n_floes; ++i )
    {
        auto const& solid = list_so.geo[i];
        auto const& movement = list_so.mov[i];

        // Import boundary
        TGeometry* geometry = new TGeometry();
        auto& boundary = geometry->outer();
        for ( auto const& point : solid.pt_own_mark.S.border )
            boundary.push_back( { point[0], point[1] } );

        // Import mesh
        TMesh* mesh = new TMesh();
        auto& points = mesh->points();
        for ( auto const& point : solid.pt_own_mark.S.vertices )
            points.push_back( { point[0], point[1] } );

        for ( auto const& triangle : solid.pt_own_mark.S.triangles )
        {
            mesh->add_triangle(
                static_cast<std::size_t>(triangle[0]) - 1,
                static_cast<std::size_t>(triangle[1]) - 1,
                static_cast<std::size_t>(triangle[2]) - 1
            );
        }

        // Create static floe
        TStaticFloe* static_floe = new TStaticFloe();
        static_floe->attach_geometry_ptr(geometry);
        static_floe->attach_mesh_ptr(mesh);
        static_floe->set_density( movement.mass_tot / static_floe->area() );

        // Import space-time state
        TState state;
        state.pos = { movement.center[0], movement.center[1] };
        state.speed = { movement.vcenter[0], movement.vcenter[1] };
        state.theta = movement.theta;
        state.rot = movement.vtheta;

        // Create Kinematic floe
        TKinematicFloe* floe = new TKinematicFloe();
        floe->attach_static_floe_ptr(static_floe);
        floe->set_state( state );
        floe->is_obstacle() = movement.obs > 0;

        // Done.
        floes.push_back(floe);
    }

    return floes;
}

}}} // namespace floe::io::matlab

#endif // FLOE_IO_MATLAB_LIST_SO_TO_FLOES_HPP
