/*!
 * \file floe/collision/matlab/periodic_detector.h
 * \brief Collision detector (matlab version) and associated functions, taking space topology (periodicity) into account
 * \see MatlabDetector for more explanations.
 * \author Quentin Jouet
 */

#ifndef FLOE_COLLISION_MATLAB_PERIODIC_DETECTOR_HPP
#define FLOE_COLLISION_MATLAB_PERIODIC_DETECTOR_HPP

#include "floe/collision/matlab/detector.h"
#include "floe/collision/periodic_contact_point.hpp"
#include "floe/collision/matlab/periodic_proximity_data.hpp"


namespace floe { namespace collision { namespace matlab
{


template <
    typename TFloeGroup,
    typename TSpaceTopology,
    typename TProximityData = PeriodicProximityData<TFloeGroup, OptimizedFloe<typename TFloeGroup::floe_type>>,
    typename TContact = PeriodicContactPoint<typename TFloeGroup::floe_type, typename TProximityData::ghost_floe_type> 
>
class PeriodicMatlabDetector : public MatlabDetector<TFloeGroup, TProximityData, TContact>
{
public:
    using base_class = MatlabDetector<TFloeGroup, TProximityData, TContact>;

    using proximity_data_type = TProximityData;
    using contact_type = TContact;

    using floe_group_type = TFloeGroup;
    using floe_type = typename floe_group_type::floe_type;
    using value_type = typename floe_type::value_type;
    using point_type = typename floe_type::point_type;
    using topology_type = TSpaceTopology;

    //! Default constructor
    PeriodicMatlabDetector() : base_class(), m_topology{nullptr} {}

    // virtual void push_back( floe_type * floe_ptr ) override
    // {
    //     base_class::push_back(floe_ptr);
    //     std::size_t floe_id = base_class::get_nb_floes() - 1;
    //     for (auto& translation : m_topology->ghosts_0())
    //     {   
    //         base_class::m_prox_data.add_ghost(floe_id, translation);
    //     }
    // }

    virtual void set_floe_group(floe_group_type const& floe_group) override {
        base_class::set_floe_group(floe_group);
        for (std::size_t floe_id = 0; floe_id<floe_group.get_floes().size(); ++floe_id){
            for (auto& translation : m_topology->ghosts_0())
            {   
                base_class::m_prox_data.add_ghost(floe_id, translation);
            }
        }
    }

    inline void set_topology(topology_type const& t) { m_topology = &t; }

private:
    topology_type const* m_topology; //!< Space topology

    void detect() override; // initialization

    virtual contact_type create_contact(std::size_t n1, std::size_t n2, point_type point1, point_type point2) const override;
};


}}} // namespace floe::collision::matlab


#endif // FLOE_COLLISION_MATLAB_PERIODIC_DETECTOR_HPP