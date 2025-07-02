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
    using real_type = typename floe_type::real_type;
    using point_type = typename floe_type::point_type;
    using topology_type = TSpaceTopology;

    //! Default constructor
    PeriodicMatlabDetector() : base_class(), m_topology{nullptr} {
        std::cout << "PeriodicMatlabDetector init" << std::endl;
    }

    virtual void set_floe_group(floe_group_type const& floe_group) override {
        base_class::set_floe_group(floe_group);
        this->add_ghosts();
    }

    virtual void rescan_floe_group() {
        base_class::rescan_floe_group();
        this->add_ghosts();
    }

    virtual void add_ghosts() {
        for (std::size_t floe_id = 0; floe_id < base_class::m_prox_data.get_floes().size(); ++floe_id){
            for (auto& translation : m_topology->ghosts_0())
            {   
                base_class::m_prox_data.add_ghost(floe_id, translation);
            }
        }
    }

    inline void set_topology(topology_type const& t) { m_topology = &t; }
    bool is_periodic() const { return true; }

private:
    topology_type const* m_topology; //!< Space topology

    void detect() override; // initialization

    virtual contact_type create_contact(std::size_t n1, std::size_t n2, point_type point1, point_type point2) const override;
};


}}} // namespace floe::collision::matlab


#endif // FLOE_COLLISION_MATLAB_PERIODIC_DETECTOR_HPP