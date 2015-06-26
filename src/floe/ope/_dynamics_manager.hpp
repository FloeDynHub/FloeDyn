/*!
 * \file ope/dynamics_manager.hpp
 * \brief Dynamics manager
 * \author Quentin Jouet
 */

#ifndef OPE_DYNAMICS_MANAGER_HPP
#define OPE_DYNAMICS_MANAGER_HPP

#include "floe/ope/FEM_manager.hpp"

namespace floe { namespace ope
{

/*! DynamicsManager
 *
 * Operator for dynamics processing (update position and speed)
 *
 */

template<typename TFloeGroup>
class DynamicsManager
{

public:
    using manager_h_type = floe::ope::FEMManager<TFloeGroup>;

    inline manager_h_type& get_manager_h(){ return m_manager_h; }

private:

    manager_h_type m_manager_h;

};

}} // namespace floe::ope


#endif // OPE_DYNAMICS_MANAGER_HPP