/*!
 * \file ope/toric_topology.hpp
 * \brief Space topology and ghost creation strategy for a periodic rectangle
 * \author Quentin Jouet
 */

#ifndef TOPOLOGY_TORIC_TOPOLOGY_HPP
#define TOPOLOGY_TORIC_TOPOLOGY_HPP

#include "floe/geometry/arithmetic/point_operators.hpp"
#include <iostream> // DEBUG

namespace floe { namespace topology
{

/*! ToricTopology
 *
 * Describes a rectangle periodic domain
 *
 */


template<
    typename TPoint
>
class ToricTopology
{
public:

    using point_list = std::vector<TPoint>;
    using T = decltype(TPoint::x);

    ToricTopology() : ToricTopology(-2, 4, -2, 4) {}
    ToricTopology(T min_x, T max_x, T min_y, T max_y) :
        m_min_x{min_x}, m_max_x{max_x},
        m_min_y{min_y}, m_max_y{max_y},
        delta_x{max_x - min_x}, delta_y{max_y - min_y} {}

    point_list ghosts(const TPoint& pt) const
    {
        return {
            pt + TPoint{delta_x, -delta_y},
            pt + TPoint{delta_x, 0},
            pt + TPoint{delta_x, delta_y},
            pt + TPoint{0, delta_y},
        };
    }

    void replace(TPoint& pt, TPoint& trans) const
    {
        while (pt.x <= m_min_x) { pt.x += delta_x; trans.x -= delta_x; }
        while (pt.x >= m_max_x) { pt.x -= delta_x; trans.x += delta_x; }
        while (pt.y <= m_min_y) { pt.y += delta_y; trans.y -= delta_y; }
        while (pt.y >= m_max_y) { pt.y -= delta_y; trans.y += delta_y; }
    }

    point_list ghosts_0() const
    {
        return {
            TPoint{delta_x, -delta_y},
            TPoint{delta_x, 0},
            TPoint{delta_x, delta_y},
            TPoint{0, delta_y},
        };
    }

private:
    T m_min_x;
    T m_max_x;
    T m_min_y;
    T m_max_y;

    T delta_x;
    T delta_y;
};


}} // namespace floe::topology


#endif // TOPOLOGY_TORIC_TOPOLOGY_HPP