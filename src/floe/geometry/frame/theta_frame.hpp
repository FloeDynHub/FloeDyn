/*!
 * \file floe/geometry/frame/theta_frame.hpp
 * \brief Frame defined by origin and the orientation of his basis
 * \author Roland Denis
 */

#ifndef FLOE_GEOMETRY_FRAME_THETA_FRAME_HPP
#define FLOE_GEOMETRY_FRAME_THETA_FRAME_HPP

#include <boost/concept/assert.hpp>
#include <cmath>

#include "floe/geometry/geometries/concepts/point_concept.hpp"
#include "floe/geometry/core/coordinate_type.hpp"

namespace floe { namespace geometry { namespace frame
{

/*! A frame of the 2D vector space, with a direct orthonormal basis
 *
 * This frame is defined by his center and angle.
 *
 * \tparam TPoint Point type
 * \tparam T      Orientation representation
 */
template <
    typename TPoint,
    typename T = typename floe::geometry::coordinate_type<TPoint>::type
>
class ThetaFrame
{
    
    BOOST_CONCEPT_ASSERT( (boost::geometry::concept::Point<TPoint>) );

public:

    // Type traits
    typedef TPoint point_type;
    typedef T      theta_type;
    typedef typename floe::geometry::coordinate_type<TPoint>::type coordinate_type;

    //! Default constructor
    ThetaFrame() : m_center{0,0}, m_theta{0} {}
    
    //! Constructor from center point and orientation
    ThetaFrame(point_type const& center, theta_type const& theta)
        : m_center{center}, m_theta{theta} {}
    
    //! Constructor from center's coordinates and orientation
    ThetaFrame(coordinate_type const& x, coordinate_type const& y, theta_type const& theta)
        : m_center{x, y}, m_theta{theta} {}

    ThetaFrame(ThetaFrame<TPoint, T> const&) = default;
    ThetaFrame& operator= (const ThetaFrame&) = default;

    //! Manipulate center
    inline TPoint& center() { return m_center; }
    inline TPoint const& center() const { return m_center; }

    //! Manipulate orientation
    inline T& theta() { return m_theta; }
    inline T const& theta() const { return m_theta; }

    //! Get first basis vector
    inline TPoint u() const
    {
        return {
            std::cos(m_theta),
            std::sin(m_theta)
        };
    }

    //! Get second basis vector
    inline TPoint v() const
    {
        return {
            - std::sin(m_theta),
            std::cos(m_theta)
        };
    }

private:
    TPoint m_center; //! center of the frame
    T m_theta; //! orientation of the frame
    
};

}}} // namespace floe::geometry::frame

#endif // FLOE_GEOMETRY_FRAME_THETA_FRAME_HPP

