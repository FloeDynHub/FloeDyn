/*!
 * \file uv_frame.hpp
 * \brief Frame defined by origin and the two basis vectors.
 * \author Roland Denis
 */

#ifndef FLOE_GEOMETRY_FRAME_UV_FRAME_HPP
#define FLOE_GEOMETRY_FRAME_UV_FRAME_HPP

#include <boost/concept/assert.hpp>
#include <cmath>

#include "floe/geometry/geometries/concepts/point_concept.hpp"
#include "floe/geometry/core/coordinate_type.hpp"
#include "floe/geometry/core/access.hpp"

namespace floe { namespace geometry { namespace frame
{

/*! A frame of the 2D vector space, with a direct orthonormal basis
 *
 * This frame is defined by his center and the first basis vector
 *
 * \tparam TPoint Point type
 * \tparam T      Orientation representation
 */
template <
    typename TPoint,
    typename T = typename floe::geometry::coordinate_type<TPoint>::type
>
class UVFrame
{
    
    BOOST_CONCEPT_ASSERT( (boost::geometry::concept::Point<TPoint>) );

public:

    // Type traits
    typedef TPoint point_type;
    typedef T      theta_type;
    typedef typename floe::geometry::coordinate_type<TPoint>::type coordinate_type;

    //! Default constructor
    UVFrame() : m_center{0,0}, m_u{1,0} {}
    
    /*! Constructor from center point and first basis vector
     * \param center    Origin of the frame.
     * \param u         First vector of the direct orthonormal basis.
     */
    UVFrame(point_type const& center, point_type const& u)
        : m_center{center}, m_u{u} {}
   
    // Center accessors
    inline TPoint& center() { return m_center; } //! Mutable center accessor.
    inline TPoint const& center() const { return m_center; } //! Constant center accessor.

    //! Get orientation
    inline T theta() const 
    { 
        using floe::geometry::get;
        return std::atan2( get<1>(m_u), get<0>(m_u) );
    }

    // First basis vector accessors
    inline TPoint const&    u() const   { return m_u; } //! Constant first basis vector accessor
    inline TPoint &         u()         { return m_u; } //! Mutable first basis vector accessor

    /*! Modify first basis vector
     * \param u The first basis vector.
     */
    inline void set_u( TPoint const& u )
    {
        using floe::geometry::get;
        using floe::geometry::set;
        set<0>(m_u, get<0>(u));
        set<1>(m_u, get<1>(u));
    }

    //! Constance second basis vector accessors
    inline TPoint   v() const   
    { 
        using floe::geometry::get;
        return { -get<1>(m_u), get<0>(m_u) };
    }

    /*! Modify second basis vector.
     * In fact, it modifies the first basis vector accordingly.
     * \param v The second basis vector.
     */
    inline void set_v( TPoint const& v )
    {
        using floe::geometry::get;
        using floe::geometry::set;
        set<0>(m_u,  get<1>(v));
        set<1>(m_u, -get<0>(v));
    }

private:
    TPoint m_center; //! center of the frame
    TPoint m_u; //! first basis vector
    
};

}}} // namespace floe::geometry::frame

#endif // FLOE_GEOMETRY_FRAME_UV_FRAME_HPP

