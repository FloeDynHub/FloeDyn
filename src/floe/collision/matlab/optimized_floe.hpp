/*!
 * \file optimized_floe.hpp
 * \brief Optimization for distance calculation and collision detection.
 * \author Roland Denis
 */

#ifndef FLOE_COLLISION_MATLAB_OPTIMIZED_FLOE_HPP
#define FLOE_COLLISION_MATLAB_OPTIMIZED_FLOE_HPP

#include <cstddef>
#include <limits>
#include <iostream> // DEBUG

// Geometry
#include "floe/geometry/geometry.hpp"
#include "floe/geometry/iterators/closing_iterator.hpp"
#include "floe/geometry/arithmetic/point_operators.hpp"
#include "floe/geometry/geometries/multi_point.hpp"

// Floe
#include "floe/floes/floe_exception.hpp"

// Frames
#include "floe/geometry/frame/uv_frame.hpp"
#include "floe/geometry/frame/frame_transformers.hpp"

// Circles
#include "floe/geometry/geometries/circle.hpp"
#include "floe/geometry/geometries/multi_circle.hpp"

#include "floe/collision/matlab/optim_interface.hpp"

namespace floe { namespace collision { namespace matlab
{

/*! Detector optimization for a floe
 *
 * It stores optimization data to accelerate distance calculation
 * and collision detection for a floe.
 *
 * \tparam TFloe Type of floe.
 */
template <
    typename TFloe
>
class OptimizedFloe : public OptimInterface<TFloe>
{

public:
    // Type traits
    using floe_type = TFloe;
    using optim_interface_type = OptimInterface<TFloe>;
    using frame_type = typename floe_type::frame_type;
    using point_type = typename optim_interface_type::point_type;
    using circle_type = typename optim_interface_type::circle_type;
    using multi_circle_type = typename optim_interface_type::multi_circle_type;
    using local_points_type = typename optim_interface_type::local_points_type;
    using real_type = typename optim_interface_type::real_type;

    /*! Constructor
     *
     * It attachs this instance to a floe and initialize optimization datas.
     *
     * \param floe Optimized floe.
     */
    OptimizedFloe( floe_type const& floe ) 
        : m_floe(floe), m_frame(floe.frame())
    {
        init();
    }

    //! Deleted default constructor
    OptimizedFloe() = delete;

    //! Default destructor
    virtual ~OptimizedFloe() = default;
    
    /*! Update optimizer
     *
     * Detects frame change of the floe and update optimization datas.
     */
    void update(bool update_local_disks=true);

    //! Global disk accessors
    inline circle_type const&   global_disk()   const   { return m_global_disk; }
    inline circle_type &        global_disk()           { return m_global_disk; }

    //! Local disks accessors
    inline multi_circle_type const& local_disks()   const   { return m_local_disks; }
    inline multi_circle_type &      local_disks()           { return m_local_disks; }

    //! Local points accessors
    inline local_points_type const& local_points()  const   { return m_local_points; }
    inline local_points_type &      local_points()          { return m_local_points; }

    real_type const&         cdist() const { return m_cdist; }
    real_type const&         tau() const { return m_tau; }
    real_type          m_cdist;        //!< Collision distance
    real_type          m_tau;          //!< Distance between surrounding disk and the floe border.

private:
    floe_type const&    m_floe; //!< Floe
    frame_type          m_frame; //!< Current frame
    circle_type         m_global_disk;  //!< The surrounding disk
    multi_circle_type   m_local_disks;  //!< The local disks (surrounding the border)
    local_points_type   m_local_points; //!< Index of first point of the border that is in the corresponding local disk

    /*! Optimizer initialization
     *
     * It creates the locals disks and the global surrounding disk.
     */
    void init();

};

/*! Return circle envelope of a geometry
 * \TODO: put that elsewhere !
 */
template < 
    typename TCircle,
    typename TMultiPoint
>
TCircle circle_envelope ( TMultiPoint const& points )
{
    using namespace floe::geometry;
    typedef typename point_type<TMultiPoint>::type point_type;
    typedef typename coordinate_type<TCircle>::type real_type;
    
    // Circle center
    const auto center = return_centroid<point_type>(points);

    // Circle radius
    point_type max_p;
    real_type max_d = 0;
    for ( const auto& pt : points )
    {
        const real_type d = comparable_distance(pt, center);
        if (d >= max_d)
        {
            max_d = d;
            max_p = pt;
        }
    }
    const real_type radius = distance(max_p, center);

    return {center, radius};
}


/*! Initializing local disks and surround disk
 */
template< typename TFloe >
void
OptimizedFloe<TFloe>::init()
{
    if (! m_floe.has_geometry() ) { throw floe::floes::FloeException("Collision optimization needs a floe with a geometry."); }
    using namespace floe::geometry;

    // Collision distance
    m_cdist = std::sqrt(m_floe.area()) / 100;

    // Get boundary
    const auto boundary = exterior_ring(m_floe.geometry());
    const std::size_t n_points = num_points(boundary);

    // Points count per disk
    const std::size_t n_pts_disk = (n_points < 50 ) ? 2 : ( (n_points < 120) ? 4 : 6 );
        
    // Initialization
    m_local_disks.resize(0);
    m_local_points.resize(0);

    // Construction of the disks
    MultiPoint<point_type> points;
    m_local_points.push_back(0);
    std::size_t cnt = 0;
    auto it = closing_iterator<decltype(boundary)>{boundary};
    const auto it_end = closing_iterator<decltype(boundary)>{boundary, true};
    for ( ; it != it_end ; ++it, ++cnt)
    {
        points.push_back(*it);

        if ( points.size() >= n_pts_disk )
        {
            m_local_disks.push_back( return_buffer<circle_type>(circle_envelope<circle_type>(points), m_cdist) );
            m_local_points.push_back( cnt );
            points.resize(0);
            points.push_back(*it);
        }
    }

    // Last local disk
    if ( points.size() > 1 )
    {
        m_local_disks.push_back( return_buffer<circle_type>(circle_envelope<circle_type>(points), m_cdist) );
        m_local_points.push_back( cnt );
    }

    //// Surrounding disk ////
    points.resize(0);
    real_type max_radius = 0.;
    for (auto const& disk : m_local_disks)
    {
        //points.push_back( disk.center );
        max_radius = std::max( max_radius, disk.radius );
    }
    m_tau = max_radius + max_radius/5; // ??? (see create_disk.m, l.16)
    //m_global_disk = return_buffer<circle_type>( circle_envelope<circle_type>(points), m_tau );
    m_global_disk = return_buffer<circle_type>( circle_envelope<circle_type>(boundary), m_tau );
}

//! Update optimization datas.
template< typename TFloe >
void
OptimizedFloe<TFloe>::update(bool update_local_disks)
{
    // New frame
    const auto new_frame = m_floe.frame();

    // Transformation from current frame to new frame
    const auto trans = geometry::frame::itransformer( m_frame, new_frame );

    // Transforming global disk
    geometry::transform( m_global_disk, m_global_disk, trans );

    if (update_local_disks){
        // Transforming local disks
        geometry::transform( m_local_disks, m_local_disks, trans );
    }

    m_frame = new_frame;
}


}}} // namespace floe::collision::matlab

#endif // FLOE_COLLISION_MATLAB_OPTIMIZED_FLOE_HPP

