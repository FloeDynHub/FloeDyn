/*!
 * \file ghost_optimized_floe.hpp
 * \brief translated reflection of an optimized floe
 * \author Quentin Jouet
 */

#ifndef FLOE_COLLISION_MATLAB_GHOST_OPTIMIZED_FLOE_HPP
#define FLOE_COLLISION_MATLAB_GHOST_OPTIMIZED_FLOE_HPP

#include <iostream> // DEBUG

#include "floe/geometry/geometry.hpp"
#include "floe/geometry/arithmetic/point_operators.hpp"
#include "floe/collision/matlab/optim_interface.hpp"

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

/*! Detector optimization for a ghost floe
 *  Refers to the real floe's optimization
 *
 * \tparam TOptim Type of Optimized floe.
 */
template <
    typename TOptim
>
class GhostOptimizedFloe : public OptimInterface<typename TOptim::floe_type>
{

public:
    // Type traits
    using optim_type = TOptim;
    using point_type = typename optim_type::point_type;
    using circle_type = typename optim_type::circle_type;
    using multi_circle_type = typename optim_type::multi_circle_type;
    using local_points_type = typename optim_type::local_points_type;
    using real_type = typename optim_type::real_type;
    using translate_strategy_type = boost::geometry::strategy::transform::translate_transformer<real_type, 2,2>;


    /*! Constructor
     * \param optim Optimized floe.
     */
    GhostOptimizedFloe( optim_type const& optim, point_type translation ) :
        m_original_id{0},
        m_optim(&optim),
        m_translation{translation},
        m_translator{translate_strategy_type(m_translation.x, m_translation.y)}
        {};

    //! Deleted default constructor
    GhostOptimizedFloe() = delete;

    //! Global disk accessor
    circle_type const& global_disk() const;
    //! Local disks accessor
    multi_circle_type const& local_disks()   const;
    //! Local points accessor
    inline local_points_type const& local_points() const { return m_optim->local_points(); }

    real_type const&         cdist() const { return m_optim->cdist(); }
    real_type const&         tau() const { return m_optim->tau(); }

    const std::size_t m_original_id; //!< Original object id in group

private:
    const optim_type* m_optim;  //!< Original object
    point_type m_translation; //!< Translation compared to original object
    translate_strategy_type m_translator; //! Translation strategy for geometry transformation

    mutable circle_type         m_global_disk;  //!< The surrounding disk
    mutable multi_circle_type   m_local_disks;  //!< The local disks (surrounding the border)
};


template<typename TOptim>
typename GhostOptimizedFloe<TOptim>::circle_type const&
GhostOptimizedFloe<TOptim>::global_disk() const 
{   
    boost::geometry::transform(m_optim->global_disk(), m_global_disk, m_translator);
    return m_global_disk;
}

template<typename TOptim>
typename GhostOptimizedFloe<TOptim>::multi_circle_type const&
GhostOptimizedFloe<TOptim>::local_disks() const 
{ 
    boost::geometry::transform(m_optim->local_disks(), m_local_disks, m_translator);
    return m_local_disks;
}


}}} // namespace floe::collision::matlab

#endif // FLOE_COLLISION_MATLAB_GHOST_OPTIMIZED_FLOE_HPP

