/*!
 * \file floes/false_hdf5_manager.hpp
 * \brief Empty HDF5 manager for g++ testing (no hdf5 support for MacOs/gcc)
 * \author Quentin Jouet
 */

#ifndef FLOE_IO_HDF5_MANAGER_HPP
#define FLOE_IO_HDF5_MANAGER_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include "floe/floes/floe_group.hpp"

namespace floe { namespace io
{

/*! HDF5Manager
 *
 * Empty version (for g++ testing)
 *
 */


template <
    typename TFloeGroup, typename TDynamicsMgr
>
class HDF5Manager
{

public:
    using floe_group_type = TFloeGroup;
    using value_type = typename TFloeGroup::value_type;
    using dynamics_mgr_type = TDynamicsMgr;

    void save_step(value_type time, const floe_group_type&, const dynamics_mgr_type&){}

    double recover_states(std::string filename, value_type time, floe_group_type&, dynamics_mgr_type&){ return 0; }

};


}} // namespace floe::io


#endif // FLOE_IO_HDF5_MANAGER_HPP
