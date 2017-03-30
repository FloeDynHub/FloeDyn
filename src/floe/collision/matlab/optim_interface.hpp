/*!
 * \file optim_interface.hpp
 * \brief Interface of an Optimized Floe
 * \author Quentin Jouet
 */

#ifndef FLOE_COLLISION_MATLAB_OPTIM_INTERFACE_HPP
#define FLOE_COLLISION_MATLAB_OPTIM_INTERFACE_HPP

 // Geometry
#include "floe/geometry/geometry.hpp"
#include "floe/geometry/iterators/closing_iterator.hpp"
#include "floe/geometry/geometries/multi_point.hpp"

// Floe
#include "floe/floes/floe_exception.hpp"

// Frames
#include "floe/geometry/frame/uv_frame.hpp"
#include "floe/geometry/frame/frame_transformers.hpp"

// Circles
#include "floe/geometry/geometries/circle.hpp"
#include "floe/geometry/geometries/multi_circle.hpp"

namespace floe { namespace collision { namespace matlab
{

/*! Minimal interface for optimized floe
 *
 *
 * \tparam TFloe Type of floe.
 */
template <
    typename TFloe
>
class OptimInterface
{

public:
    // Type traits
    using floe_type = TFloe;
    using point_type = typename floe_type::point_type;
    using circle_type = floe::geometry::Circle<point_type>;
    using multi_circle_type = floe::geometry::MultiCircle<circle_type>;
    using local_points_type = std::vector<std::size_t>;
    using real_type = typename floe_type::real_type;

    //! Global disk accessor
    virtual circle_type const& global_disk() const = 0;
    //! Local disks accessor
    virtual multi_circle_type const& local_disks()   const = 0;
    //! Local points accessor
    virtual local_points_type const& local_points() const = 0;

    //! Contact distance accessor
    virtual real_type const&         cdist() const = 0;
    //! Accessor for distance between surrounding disk and the floe border.
    virtual real_type const&         tau() const = 0;

};


}}} // namespace floe::collision::matlab

#endif // FLOE_COLLISION_MATLAB_OPTIM  _INTERFACE_HPP

