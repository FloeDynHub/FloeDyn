/*!
 * \file floe/collision/matlab/periodic_detector.h
 * \brief Collision detector (matlab version) and associated functions, taking space topology (periodicity) into account
 * \see MatlabDetector for more explanations.
 * \author Quentin Jouet
 */

#ifndef FLOE_COLLISION_MATLAB_PERIODIC_DETECTOR_HPP
#define FLOE_COLLISION_MATLAB_PERIODIC_DETECTOR_HPP

#include "floe/collision/matlab/detector.h"
#include "floe/collision/matlab/ghost_optimized_floe.hpp"
#include "floe/floes/ghost_floe.hpp"
#include "floe/collision/periodic_contact_point.hpp"


namespace floe { namespace collision { namespace matlab
{


template <
    typename TFloe,
    typename TSpaceTopology,
    typename TGhostFloe = floes::GhostFloe<TFloe>,
    typename TContact = PeriodicContactPoint<TFloe, TGhostFloe> 
>
class PeriodicMatlabDetector : public MatlabDetector<TFloe, TContact>
{
public:
    using base_class = MatlabDetector<TFloe, TContact>;

    using floe_type = typename base_class::floe_type;
    using optim_type = typename base_class::optim_type;
    using value_type = typename base_class::value_type;
    using floe_interface_type = typename base_class::floe_interface_type;
    using optim_interface_type = typename base_class::optim_interface_type;
    using point_type = typename base_class::point_type;
    using contact_type = TContact;
    using ghost_floe_type = TGhostFloe;
    using ghost_optim_type = GhostOptimizedFloe<optim_type>;
    using topology_type = TSpaceTopology;

    //! Default constructor
    PeriodicMatlabDetector() : base_class(), m_topology{nullptr} {}

    // PeriodicMatlabDetector(topology_type& topology) : base_class(), m_topology{topology} {}

    virtual void push_back( floe_type * floe_ptr )
    {
        base_class::push_back(floe_ptr);
        std::size_t floe_id = base_class::m_floes.size() - 1;
        optim_type* optim_ptr = base_class::m_optims[floe_id];
        for (auto& translation : m_topology->ghosts_0())
        {
            m_ghost_floes.push_back(ghost_floe_type{ *floe_ptr, translation, floe_id});
            m_ghost_optims.push_back(ghost_optim_type{ *optim_ptr, translation });
        }
    }

    //!Empty floe and optim lists
    virtual void reset() override {
        base_class::reset();
        m_ghost_floes.clear(); m_ghost_optims.clear();
    }

    virtual floe_interface_type const& get_floe(std::size_t n) const;
    virtual optim_interface_type const& get_optim(std::size_t n) const;

    inline void set_topology(topology_type const& t) { m_topology = &t; }

private:
    std::vector<ghost_floe_type>     m_ghost_floes; //!< Ghost floes list.
    std::vector<ghost_optim_type>    m_ghost_optims; //!< Ghost Optimized floes list.
    topology_type const* m_topology; //!< Space topology

    void detect(); // initialization

    inline virtual void set_dist_secu(std::size_t n1, std::size_t n2, value_type val) override {
        (n2 >= base_class::m_floes.size()) ? base_class::m_dist_secu(n1, n2) = val : base_class::m_dist_secu(n1, n2) = base_class::m_dist_secu(n2, n1) = val;
    }
    inline virtual void set_indic(std::size_t n1, std::size_t n2, short val) override {
        (n2 >= base_class::m_floes.size()) ? base_class::m_indic(n1, n2) = val : base_class::m_indic(n1, n2) = base_class::m_indic(n2, n1) = val;
    }
    virtual void set_dist_opt(std::size_t n1, std::size_t n2, value_type val);
    virtual contact_type create_contact(std::size_t n1, std::size_t n2, point_type point1, point_type point2) const;
    //! returns real floe id if param id refers to a ghost
    virtual std::size_t real_floe_id(std::size_t n) const {
        std::size_t N { base_class::m_floes.size() };
        if (n >= N)
            return m_ghost_floes[n - N].m_original_id;
        else
            return n;
    }

    friend class ope::TimeScaleManager<PeriodicMatlabDetector<TFloe, TSpaceTopology, TGhostFloe, TContact>>;
};


}}} // namespace floe::collision::matlab


#endif // FLOE_COLLISION_MATLAB_PERIODIC_DETECTOR_HPP