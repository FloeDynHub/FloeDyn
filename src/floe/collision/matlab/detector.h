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
 */
template <
    typename TFloe,
    typename TContact = ContactPoint<TFloe> 
>
class MatlabDetector
{

public:
    // Type traits
    using floe_type = TFloe;
    using optim_type = OptimizedFloe<TFloe>;
    using point_type = typename floe_type::point_type;
    using value_type = typename floe_type::value_type;
    using floe_interface_type = typename floe_type::floe_interface_type;
    using optim_interface_type = typename optim_type::optim_interface_type;
    using circle_type = typename optim_type::circle_type;
    using multi_circle_type = typename optim_type::multi_circle_type;
    using contact_type = TContact;
    using contact_graph_type = ContactGraph<contact_type>;
    using contact_list_type = typename contact_graph_type::edge_property_type::base_class;

    // typedef ublas::symmetric_matrix<value_type, ublas::lower> dist_matrix_type; // Type of distance matrix
    // typedef ublas::symmetric_matrix<std::size_t, ublas::lower> indic_matrix_type; // Type of indicator matrix
    using dist_matrix_type = ublas::matrix<value_type>; //!< Type of distance matrix
    using indic_matrix_type = ublas::matrix<short>; //!< Type of indicator matrix

    //! Default constructor
    MatlabDetector()
        : m_floes{}, m_optims{}, 
          m_indic{0,0}, m_dist_secu{0,0}, m_dist_opt{0,0}, m_detection_mode{0}, m_detection_chgt{1}
    {}

    //! Deleted copy constructor
    MatlabDetector( MatlabDetector<TFloe, TContact> const& ) = delete;

    //! Deleted copy operator
    MatlabDetector<TFloe, TContact>& operator= (MatlabDetector const&) = delete;

    //! Destructor
    ~MatlabDetector() { for ( auto& optim_ptr : m_optims ) delete optim_ptr; }

    /*! Add a floe in the detector scope
     * It automatically creates the optimization datas associated to the new floe.
     * \param floe_ptr Pointer to the floe to add.
     */
    virtual void push_back( floe_type * floe_ptr )
    {
        m_floes.push_back(floe_ptr);
        m_optims.push_back( new optim_type{*floe_ptr} );
        m_previous_step_states.resize(m_floes.size());
    }

    //!Empty floe and optim lists
    virtual void reset() { m_floes.clear(); m_optims.clear(); }

    /*! Update collision informations
     *
     * It updates optimization datas of all associated floes and launch contact detection.
     */
    bool update();

    //! Access contacts graph
    contact_graph_type const& contact_graph() const { return m_contacts; }

    // Some informations
    std::size_t num_local_disks() const { std::size_t cnt = 0; for (auto const& opt : m_optims) cnt += opt->local_disks().size(); return cnt; }
    std::size_t num_points() const { std::size_t cnt = 0; for (auto const& floe : m_floes) cnt += floe->geometry().outer().size(); return cnt; }

    inline bool interpenetration() const { return m_interpenetration; }

    //! Container accessors
    inline virtual floe_interface_type const& get_floe(std::size_t n) const { return *(m_floes[n]); }
    inline virtual optim_interface_type const& get_optim(std::size_t n) const { return *(m_optims[n]); }

    //! if contact has not been solved, annul dist_opt
    void clean_dist_opt();

    //! is there any floe interpenetration ? returns true if not.
    bool check_interpenetration();
    void backup_step_states();
    void recover_previous_step_states();


protected:
    std::vector<floe_type *>     m_floes; //!< Floes list.
    std::vector<optim_type*>    m_optims; //!< Optimization datas list.

    indic_matrix_type m_indic; //!< Indicator of collision (0=far away, 1=close, 2=contact)
    dist_matrix_type m_dist_secu; //!< Security distance
    dist_matrix_type m_dist_opt; //!< Optimial distance
    contact_graph_type m_contacts; //!< Contact graph
    bool m_detection_mode; //! Detection mode ('eta_min' in matlab)
    bool m_detection_chgt; //! Detection status ('eta_chgt' in matlab)
    bool m_interpenetration; //! Floe interpenetration

    std::vector<typename floe_type::state_type> m_previous_step_states;

    //! Talking about segments
    //! \todo put that somewhere else !
    typedef std::pair<point_type,point_type> segment_type;
    value_type segment_pos( segment_type const& segment, point_type const& point ) const;
    value_type segment_dist( segment_type const& segment, point_type const& point ) const; 
    segment_type segment_from_id1( std::size_t n, std::size_t id1 ) const;
    segment_type segment_from_id2( std::size_t n, std::size_t id2 ) const;
    point_type point_from_id( std::size_t n, std::size_t id ) const;
    inline point_type point_from_pos( segment_type const& segment, value_type pos ) const;

    void detection_mode(); // detection mode (min or max collision distance)
    void prepare_detection(); // preparation
    virtual void detect(); // initialization + detection
    //! Detects collisions in 4 main steps
    void detect_step1();
    void detect_step2( std::size_t n1, std::size_t n2 );
    void detect_step3( std::size_t n1, std::size_t n2, std::vector<std::size_t> const& ldisks1, std::vector<std::size_t> const& ldisks2 );
    
    template <typename TAdjacency>
    value_type detect_step4( std::size_t n1, std::size_t n2, std::vector<std::size_t> const& ldisks1, std::vector<std::size_t> const& ldisks2, TAdjacency const& adjacency);

    inline virtual void set_dist_secu(std::size_t n1, std::size_t n2, value_type val) { m_dist_secu(n1, n2) = m_dist_secu(n2, n1) = val; }
    inline virtual void set_indic(std::size_t n1, std::size_t n2, short val) { m_indic(n1, n2) = m_indic(n2, n1) = val; }
    inline virtual void set_dist_opt(std::size_t n1, std::size_t n2, value_type val) { m_dist_opt(n1, n2) = m_dist_opt(n2, n1) = val; }
    inline virtual contact_type create_contact(std::size_t n1, std::size_t n2, point_type point1, point_type point2) const {
        return { m_floes[n1], m_floes[n2], point1, point2 }; }
    inline virtual std::size_t real_floe_id(std::size_t n) const { return n; }

    friend class domain::TimeScaleManager<MatlabDetector<TFloe, TContact>>;

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
