/*!
 * \file collision/collision_manager.hpp
 * \brief Collision manager
 * \author Quentin Jouet
 */

#ifndef OPE_COLLISION_MANAGER_HPP
#define OPE_COLLISION_MANAGER_HPP

namespace floe { namespace collision
{

/*! CollisionManager
 *
 * Operator for collision processing
 *
 */

template<typename TManager_h>
class CollisionManager
{

public:
	using manager_h_type = TManager_h;

    inline manager_h_type& get_manager_h(){ return m_manager_h; }

private:

    manager_h_type m_manager_h; //!< Discrete collision manager

};

}} // namespace floe::collision


#endif // OPE_COLLISION_MANAGER_HPP