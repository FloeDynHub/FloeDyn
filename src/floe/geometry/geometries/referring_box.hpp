#ifndef FLOE_GEOMETRY_GEOMETRIES_REFERRING_BOX_HPP
#define FLOE_GEOMETRY_GEOMETRIES_REFERRING_BOX_HPP

#include <boost/geometry/geometries/box.hpp>

namespace floe { namespace geometry {

/*! Axis Aligned Box type thas is linked to an object
 *
 * \tparam TPoint   Point type
 */

template < 
    typename TPoint
    //typename Geometry
>
class ReferringBox : public Box<TPoint> {
public:
    //ReferringBox() = delete;
    //ReferringBox(Geometry const& geometry) : 

private:
    //const Geometry & m_geometry;
};

}} // namespace floe::geometry

#endif // FLOE_GEOMETRY_GEOMETRIES_BOX_HPP
