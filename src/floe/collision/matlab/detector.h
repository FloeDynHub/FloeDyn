/*!
 * \file floe/collision/matlab/detector.h
 * \brief Collision detector (matlab version) and associated functions.
 * \see MatlabDetector for more explanations.
 * \author Roland Denis, Quentin Jouet
 */

#ifndef FLOE_COLLISION_MATLAB_DETECTOR_H
#define FLOE_COLLISION_MATLAB_DETECTOR_H

#include <vector>
#include <cstddef>
#include <limits>
#include <utility>

#include <iostream> // DEBUG
#include <boost/timer/timer.hpp> // DEBUG

// uBlas
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "floe/collision/matlab/optimized_floe.hpp"
#include "floe/collision/contact_point.hpp"
#include "floe/collision/contact_graph.hpp"

#include "floe/geometry/geometry.hpp"
#include "floe/geometry/arithmetic/point_operators.hpp"
#include <boost/geometry/algorithms/intersects.hpp>
#include <algorithm>

#include "floe/domain/time_scale_manager.hpp"

#include "floe/collision/matlab/proximity_data.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace floe { namespace collision { namespace matlab
{

namespace ublas = boost::numeric::ublas;

/*! Collision detector (matlab version)
 *
 * This detector behaves like the matlab code.
 * It is only a "draft version" of what a generic detector should be.
 * We can a many detectors classes that can be put together to tune the
 * detector behavior.
 * Here we first use a global disk intersection algorithm and then a more precise 
 * algorithm to find contact point between any polygonal shapes (like search_geo_contact.m)
 * but we can imagine to use a AABB algorithm to accelerate this first step before the 
 * precise step.
 *
 * There is also a lot of functions in this file that have nothing to do here.
 *
 * \tparam TFloe Type of a floe (it must be a kinematic object).
 * \tparam TProximityData  Type of internal data.
 * \tparam TContact  Type of a contact point.
 */
template <
    typename TFloeGroup,
    typename TProximityData = ProximityData<TFloeGroup, OptimizedFloe<typename TFloeGroup::floe_type>>,
    typename TContact = ContactPoint<typename TFloeGroup::floe_type> 
>
class MatlabDetector
{

public:
    // Type traits
    using floe_group_type = TFloeGroup;
    using floe_type = typename floe_group_type::floe_type;
    using optim_type = typename TProximityData::optim_type;
    using point_type = typename floe_type::point_type;
    using real_type = typename floe_type::real_type;
    using floe_interface_type = typename floe_type::floe_interface_type;
    using optim_interface_type = typename optim_type::optim_interface_type;
    using circle_type = typename optim_type::circle_type;
    using multi_circle_type = typename optim_type::multi_circle_type;
    using contact_type = TContact;
    using contact_graph_type = ContactGraph<contact_type>;
    using contact_list_type = typename contact_graph_type::edge_property_type::base_class;
    using proximity_data_type = TProximityData;

    //! Default constructor
    MatlabDetector()
        : m_prox_data{}, m_detection_mode{0}, m_detection_chgt{1} {}

    //! Deleted copy constructor
    MatlabDetector( MatlabDetector<TFloeGroup, TContact> const& ) = delete;

    //! Deleted copy operator
    MatlabDetector<TFloeGroup, TContact>& operator= (MatlabDetector const&) = delete;

    //! Destructor
    ~MatlabDetector() { for ( auto& optim_ptr : m_prox_data.get_optims() ) delete optim_ptr; }

    /*! Add a floe in the detector scope
     * It automatically creates the optimization datas associated to the new floe.
     * \param floe_ptr Pointer to the floe to add.
     */
    // virtual void push_back( floe_type * floe_ptr )
    // {
    //     m_prox_data.push_back(floe_ptr);
    // }

    virtual void set_floe_group(floe_group_type const& floe_group){
        m_prox_data.set_floe_group(floe_group);
    }

    virtual void rescan_floe_group() {
        m_prox_data.update_optim_vars();
        m_prox_data.auto_resize();
    }

    //!Empty floe and optim lists
    virtual void reset() { m_prox_data.reset(); }

    /*! Update collision informations
     *
     * It updates optimization datas of all associated floes and launch contact detection.
     */
    bool update();

    //! Access contacts graph
    contact_graph_type const& contact_graph() const { return m_contacts; }

    // Some informations
    std::size_t num_local_disks() const { std::size_t cnt = 0; for (auto const& opt : m_prox_data.get_optims()) cnt += opt->local_disks().size(); return cnt; }
    std::size_t num_points() const { std::size_t cnt = 0; for (auto const& floe : m_prox_data.get_floes() ) cnt += floe.geometry().outer().size(); return cnt; }

    //! Container accessors
    inline std::size_t get_nb_floes() const { return m_prox_data.nb_floes(); }
    inline floe_type const& get_floe(std::size_t n) const { return m_prox_data.get_floe(n); }
    inline optim_type const& get_optim(std::size_t n) const { return m_prox_data.get_optim(n); }
    inline virtual floe_interface_type const& get_floe_itf(std::size_t n) const { return m_prox_data.get_floe_itf(n); }
    inline virtual optim_interface_type const& get_optim_itf(std::size_t n) const { return m_prox_data.get_optim_itf(n); }

    //! if contact has not been solved, cancel its dist_opt
    void clean_dist_opt();

    //! is there any floe interpenetration ? returns true if not.
    bool check_interpenetration();

    inline proximity_data_type const& data() { return m_prox_data; }

protected:
    proximity_data_type m_prox_data;
    contact_graph_type m_contacts; //!< Contact graph
    bool m_detection_mode; //! Detection mode ('eta_min' in matlab)
    bool m_detection_chgt; //! Detection status ('eta_chgt' in matlab)

    inline optim_type& get_optim(std::size_t n) { return m_prox_data.get_optim(n); }

    //! Talking about segments
    //! \todo put that somewhere else !
    typedef std::pair<point_type,point_type> segment_type;
    real_type segment_pos( segment_type const& segment, point_type const& point ) const;
    real_type segment_dist( segment_type const& segment, point_type const& point ) const; 
    segment_type segment_from_id1( std::size_t n, std::size_t id1 ) const;
    segment_type segment_from_id2( std::size_t n, std::size_t id2 ) const;
    point_type point_from_id( std::size_t n, std::size_t id ) const;
    inline point_type point_from_pos( segment_type const& segment, real_type pos ) const;

    void detection_mode(); // detection mode (min or max collision distance)
    virtual void prepare_detection(); // preparation
    virtual void prepare_optims();
    void prepare_contact_graph();
    virtual void detect(); // initialization + detection
    //! Detects collisions in 4 main steps
    virtual void detect_step1();
    void detect_step2( std::size_t n1, std::size_t n2 );
    void detect_step3( std::size_t n1, std::size_t n2, std::vector<std::size_t> const& ldisks1, std::vector<std::size_t> const& ldisks2 );
    
    template <typename TAdjacency>
    real_type detect_step4( std::size_t n1, std::size_t n2, std::vector<std::size_t> const& ldisks1, std::vector<std::size_t> const& ldisks2, TAdjacency const& adjacency);

    inline virtual contact_type create_contact(std::size_t n1, std::size_t n2, point_type point1, point_type point2) const {
        return { &m_prox_data.get_floe(n1), &m_prox_data.get_floe(n2), point1, point2 }; }

};

//! Some utilities (workaround for distance to circle)
template <typename TCircle1, typename TCircle2>
inline
typename geometry::coordinate_type<TCircle1>::type
distance_circle_circle( TCircle1 const& c1, TCircle2 const& c2 )
{
    using namespace floe::geometry;
    return distance( center_view<const TCircle1>(c1), center_view<const TCircle2>(c2)) - get_radius(c1) - get_radius(c2);
}

template <typename TPoint, typename TCircle>
inline
typename geometry::coordinate_type<TPoint>::type
distance_point_circle( TPoint const& point, TCircle const& circle )
{
    using namespace floe::geometry;
    return distance(point, center_view<const TCircle>(circle)) - get_radius(circle);
}


}}} // namespace floe::collision::matlab

#endif // FLOE_COLLISION_MATLAB_DETECTOR_H
