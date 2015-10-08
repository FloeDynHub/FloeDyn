/*!
 * \file ope/collision_manager.hpp
 * \brief Collision manager
 * \author Quentin Jouet
 */

#ifndef OPE_COLLISION_MANAGER_HPP
#define OPE_COLLISION_MANAGER_HPP

#include "floe/ope/LCP_manager.hpp"

namespace floe { namespace ope
{

/*! CollisionManager
 *
 * Operator for collision processing
 *
 */

template<typename TFloe>
class CollisionManager
{

public:
    using manager_h_type = floe::ope::LCPManager<typename TFloe::value_type>;

    inline manager_h_type& get_manager_h(){ return m_manager_h; }

private:

    manager_h_type m_manager_h;

};

}} // namespace floe::ope


#endif // OPE_COLLISION_MANAGER_HPP