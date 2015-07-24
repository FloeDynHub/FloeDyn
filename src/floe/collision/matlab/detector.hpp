/*!
 * \file floe/collision/matlab/detector.hpp
 * \brief Collision detector (matlab version) and associated functions.
 * \see MatlabDetector for more explanations.
 * \author Roland Denis, Quentin Jouet
 */

#ifndef FLOE_COLLISION_MATLAB_DETECTOR_HPP
#define FLOE_COLLISION_MATLAB_DETECTOR_HPP

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
    using indic_matrix_type = ublas::matrix<std::size_t>; //!< Type of indicator matrix

    //! Default constructor
    MatlabDetector()
        : m_floes{}, m_optims{}, 
          m_indic{0,0}, m_dist_secu{0,0}, m_dist_opt{0,0}
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
    virtual void push_back( floe_type const* floe_ptr )
    {
        m_floes.push_back(floe_ptr);
        m_optims.push_back( new optim_type{*floe_ptr} );
    }

    /*! Update collision informations
     *
     * It updates optimization datas of all associated floes and launch contact detection.
     */
    void update();

    //! Access contacts graph
    contact_graph_type const& contact_graph() const { return m_contacts; }

    // Some informations
    std::size_t num_local_disks() const { std::size_t cnt = 0; for (auto const& opt : m_optims) cnt += opt->local_disks().size(); return cnt; }
    std::size_t num_points() const { std::size_t cnt = 0; for (auto const& floe : m_floes) cnt += floe->geometry().outer().size(); return cnt; }

    // Accessors for time_scale_manager
    inline dist_matrix_type const& get_dist_secu() const { return m_dist_secu; }
    inline dist_matrix_type const& get_dist_opt() const { return m_dist_opt; }
    inline indic_matrix_type const& get_indic() const { return m_indic; }
    inline std::vector<floe_type const*> const& get_floes() const { return m_floes; }
    inline std::vector<optim_type*> const& get_optims() const { return m_optims; }
    // End Accessors for time_scale_manager

    //! Container accessors
    inline virtual floe_interface_type const& get_floe(std::size_t n) const { return *(m_floes[n]); }
    inline virtual optim_interface_type const& get_optim(std::size_t n) const { return *(m_optims[n]); }

    //! if contact has not been solved, annul dist_opt
    void clean_dist_opt();

protected:
    std::vector<floe_type const*>     m_floes; //!< Floes list.
    std::vector<optim_type*>    m_optims; //!< Optimization datas list.

    indic_matrix_type m_indic; //!< Indicator of collision (0=far away, 1=close, 2=contact)
    dist_matrix_type m_dist_secu; //!< Security distance
    dist_matrix_type m_dist_opt; //!< Optimial distance
    contact_graph_type m_contacts; //!< Contact graph

    //! Talking about segments
    //! \todo put that somewhere else !
    typedef std::pair<point_type,point_type> segment_type;
    value_type segment_pos( segment_type const& segment, point_type const& point ) const;
    value_type segment_dist( segment_type const& segment, point_type const& point ) const; 
    segment_type segment_from_id1( std::size_t n, std::size_t id1 ) const;
    segment_type segment_from_id2( std::size_t n, std::size_t id2 ) const;
    point_type point_from_id( std::size_t n, std::size_t id ) const;
    inline point_type point_from_pos( segment_type const& segment, value_type pos ) const;

    void prepare_detection(); // preparation
    virtual void detect(); // initialization + detection
    //! Detects collisions in 4 main steps
    void detect_step1();
    void detect_step2( std::size_t n1, std::size_t n2 );
    void detect_step3( std::size_t n1, std::size_t n2, std::vector<std::size_t> const& ldisks1, std::vector<std::size_t> const& ldisks2 );
    
    template <typename TAdjacency>
    value_type detect_step4( std::size_t n1, std::size_t n2, std::vector<std::size_t> const& ldisks1, std::vector<std::size_t> const& ldisks2, TAdjacency const& adjacency);

    inline virtual void set_dist_secu(std::size_t n1, std::size_t n2, value_type val) { m_dist_secu(n1, n2) = m_dist_secu(n2, n1) = val; }
    inline virtual void set_indic(std::size_t n1, std::size_t n2, std::size_t val) { m_indic(n1, n2) = m_indic(n2, n1) = val; }
    inline virtual void set_dist_opt(std::size_t n1, std::size_t n2, value_type val) { m_dist_opt(n1, n2) = m_dist_opt(n2, n1) = val; }
    inline virtual contact_type create_contact(std::size_t n1, std::size_t n2, point_type point1, point_type point2) const {
        return { m_floes[n1], m_floes[n2], point1, point2 }; }
    inline virtual std::size_t real_floe_id(std::size_t n) const { return n; }

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


//! Update detector
template <
    typename TFloe,
    typename TContact
>
void
MatlabDetector<TFloe, TContact>::update()
{   
    // update data and clear contacts
    prepare_detection();
    // Launch collisions detection
    detect();
}


//! preparing members for new detection
template <
    typename TFloe,
    typename TContact
>
void
MatlabDetector<TFloe, TContact>::prepare_detection()
{
    // Update floe optimizers
    for ( auto optim_ptr : m_optims )
        optim_ptr->update();

    // Initializing graph
    m_contacts.clear();
    for ( auto floe_ptr : m_floes )
        add_vertex(floe_ptr, m_contacts);

}


//! Starting detection
template <
    typename TFloe,
    typename TContact
>
void
MatlabDetector<TFloe, TContact>::detect()
{

    // Number of floes
    const std::size_t N = m_floes.size();

    // Resize matrix
    m_indic.resize(N, N);
    // m_dist_secu.resize(N, N);
    // m_dist_opt.resize(N, N);
    m_dist_opt = ublas::scalar_matrix<value_type>(N, N, 0);
    m_dist_secu = ublas::scalar_matrix<value_type>(N, N, 0);

    std::cout << N << " " << std::endl;

    // Detection
    detect_step1();
}


//! Tests global disk intersection
template <
    typename TFloe,
    typename TContact
>
void
MatlabDetector<TFloe, TContact>::detect_step1()
{

    // Number of floes
    const std::size_t N = m_floes.size();

    // Level 1 loop
    // TODO: intersects -like for multi_circle !!

    for (std::size_t n1 = 0; n1 < N; ++n1)
    {
        auto const& opt1 = get_optim(n1);
        for (std::size_t n2 = n1 + 1; n2 < m_dist_secu.size2(); ++n2)
        {
            auto const& opt2 = get_optim(n2);
     
    /* old version for comparaison (dev) */
    // for (std::size_t n1 = 1; n1 < N; ++n1)
    // {
    //     auto const& opt1 = get_optim(n1);
    //     for (std::size_t n2 = 0; n2 < n1; ++n2)
    //     {
    //         auto const& opt2 = get_optim(n2);

            const auto dist = distance_circle_circle( 
                opt1.global_disk(),
                opt2.global_disk()
            ) + opt1.tau() + opt2.tau();

            set_dist_secu(n1, n2, dist);

            if ( dist > std::max( opt1.cdist(), opt2.cdist() ) )
            {
                set_indic(n1, n2, 0);
            } 
            else 
            {
                detect_step2(n1, n2);
            }
        }
    }
}

//! Finds local disks that are in the other floe global disk
template <
    typename TFloe,
    typename TContact
>
void
MatlabDetector<TFloe, TContact>::detect_step2( std::size_t n1, std::size_t n2 )
{
    auto const& opt1 = get_optim(n1);
    auto const& opt2 = get_optim(n2);

    set_indic(n1, n2, 1);
    const value_type dzone = -( m_dist_secu(n1, n2) - opt1.tau() - opt2.tau() );

    // Local disks of obj1 that intersects global disk of obj2
    std::vector<std::size_t> ldisks1;
    if ( dzone > opt1.tau() )
    {
        for ( std::size_t i = 0; i < opt1.local_disks().size(); ++i )
        {
            if ( distance_circle_circle( opt1.local_disks()[i], opt2.global_disk() ) < 0 )
                ldisks1.push_back(i);
        }
    }

    // Local disks of obj2 that intersects global disk of obj1
    std::vector<std::size_t> ldisks2;
    if ( dzone > opt2.tau() )
    {
        for ( std::size_t i = 0; i < opt2.local_disks().size(); ++i )
        {
            if ( distance_circle_circle( opt2.local_disks()[i], opt1.global_disk() ) < 0 )
                ldisks2.push_back(i);
        }
    }

    // What's up doctor ?
    if ( ldisks1.size() == 0 && ldisks2.size() == 0 )
        set_dist_secu(n1, n2, dzone);
    else if ( ldisks1.size() != 0 && ldisks2.size() == 0 )
        set_dist_secu(n1, n2, opt1.tau());
    else if ( ldisks1.size() == 0 && ldisks2.size() != 0 )
        set_dist_secu(n1, n2, opt2.tau());
    else
        detect_step3(n1, n2, ldisks1, ldisks2);
}

//! Adjacency matrix of local disks
template <
    typename TFloe,
    typename TContact
>
void
MatlabDetector<TFloe, TContact>::detect_step3( 
    std::size_t n1, std::size_t n2, 
    std::vector<std::size_t> const& ldisks1, std::vector<std::size_t> const& ldisks2 
)
{
    auto const& opt1 = get_optim(n1);
    auto const& opt2 = get_optim(n2);
    
    // Adjacency matrix
    boost::numeric::ublas::compressed_matrix<bool> adjacency(ldisks1.size(), ldisks2.size());

    // Intersection counter
    std::size_t cnt = 0;

    // Security and optimal distance
    value_type dist_s = std::numeric_limits<value_type>::max();
    value_type dist_o = dist_s;


    // Intersection finding
    for ( std::size_t i = 0; i < ldisks1.size(); ++i )
    {
        const circle_type d1 = opt1.local_disks()[ldisks1[i]];
        for ( std::size_t j = 0; j < ldisks2.size(); ++j )
        {
            const circle_type d2 = opt2.local_disks()[ldisks2[j]];
            const value_type dist = distance_circle_circle( d1, d2 );
            dist_o = std::min( dist_o, dist + std::min(d1.radius, d2.radius) ); // Not exactly like in the matlab code // quentin : can be optimized

            if ( dist < 0 )
            {
                adjacency(i,j) = 1;
                dist_s = std::min( dist_s, dist + opt1.tau() + opt2.tau() );
                ++cnt;
            }
        }
    }

    // What's up doctor ?
    if (cnt == 0)
    {
        set_dist_secu(n1, n2, dist_s);
        set_dist_opt(n1, n2, dist_o);
    } 
    else 
    {
        set_dist_secu(n1, n2,  std::min(
            detect_step4( n2, n1, ldisks2, ldisks1, ublas::trans(adjacency) ), // Detects contacts obj2 -> obj1
            detect_step4( n1, n2, ldisks1, ldisks2, adjacency ) // Detects contacts obj1 -> obj2
        ));
    }
}

//! Searching contacts
template <
    typename TFloe,
    typename TContact
>
template <typename TAdjacency>
typename MatlabDetector<TFloe, TContact>::value_type
MatlabDetector<TFloe, TContact>::
detect_step4( 
    std::size_t n1, std::size_t n2, 
    std::vector<std::size_t> const& ldisks1, std::vector<std::size_t> const& ldisks2,
    TAdjacency const& adjacency
)
{
    using namespace floe::geometry;

    auto const& opt1 = get_optim(n1);
    auto const& opt2 = get_optim(n2);

    // Contact list
    contact_list_type contact_list;

    value_type global_min_dist = std::numeric_limits<value_type>::max(); // Minimum distance from any points of obj1 to obj2

    // !!! TEST !!!
    // std::vector<unsigned char> segment_cnt(m_floes[n2]->static_floe().geometry().outer().size(), 0);
    
    // Loop over disks of obj1
    for ( auto it1 = adjacency.begin1(); it1 != adjacency.end1(); ++it1 )
    {
        const std::size_t id1 = ldisks1[it1.index1()];

        // Loop over points of this disk
        for ( std::size_t ipt1 = opt1.local_points()[id1]; ipt1 < opt1.local_points()[id1+1]; ++ipt1 )
        {
            const point_type point1 = point_from_id(n1, ipt1);

            // Best contact
            value_type min_dist = std::numeric_limits<value_type>::max(); // Minimum distance from this point to the other floe 
            contact_type min_contact;
            // int seg_id = -1; // !!!! TEST !!!!

            bool dangling_point = false;    // Indicate a point that is mayby in the sub-derivative of the other floe
            std::size_t dangling_id = 0;    // Id of the point around which there is a dangling point
            std::size_t last_id2 = 0;       // Last disk id, to detect discontinuity in the disk list

            // Loop over disks of obj2 that intersects the disk of obj1
            for ( auto it2 = it1.begin(); it2 != it1.end(); ++it2 )
            {
                const std::size_t id2 = ldisks2[it2.index2()];
                
                if (id2 != last_id2+1) dangling_point = false; // Discontinuity in the disk list
                last_id2 = id2;
                
                // If the point is in this disk, there may be a contact
                if ( distance_point_circle(point1, opt2.local_disks()[id2]) <= 0 )
                {
                    // Loop over points of this disk
                    for ( std::size_t ipt2 = opt2.local_points()[id2]; ipt2 < opt2.local_points()[id2+1]; ++ipt2 )
                    {
                        // Segment beginning with ipt2
                        const segment_type segment = segment_from_id1(n2, ipt2);
                        
                        // Position of the projection of the point of obj1 on the segment of obj2
                        const value_type pos = segment_pos(segment, point1);
 
                        if (pos <= 0) // Backward position, handle only if the point is dangling
                        {
                            if (dangling_point)
                            {
                                // Contact point-point (check sub-derivative)
                                const point_type point2 = point_from_id(n2, ipt2);
                                const value_type dist = distance(point1, point2);
                                if (
                                        dist < min_dist
                                    &&  segment_pos( segment_from_id2(n1, ipt1), point2 ) >= 1
                                    &&  segment_pos( segment_from_id1(n1, ipt1), point2 ) <= 0
                                )
                                {
                                    min_contact = create_contact(n1, n2, point1, point2);
                                    min_dist = dist;
                                    // seg_id = -1; // TEST
                                }

                                dangling_point = false;
                            }
                        }
                        else if (pos > 0 && pos < 1) // Middle position => point-segment contact
                        {
                            // Contact point-segment
                            const point_type point2 = point_from_pos(segment, pos);
                            const value_type dist = distance(point1, point2);
                            if ( dist < min_dist )
                            {
                                min_contact = create_contact(n1, n2, point1, point2);
                                // min_contact = { m_floes[n2], m_floes[n1], point2, point1 }; // TEST
                                min_dist = dist;
                                // seg_id = ipt2; // TEST
                            }
                            dangling_point = false; 
                        } 
                        else if (pos >= 1) // Forward position, mark the point as dangling
                        {
                            dangling_point = true;
                            dangling_id = ipt2+1;
                        }

                    } // Loop over points of this disk
                } 
                else
                {   
                    // TODO : synchro matlab
                    const value_type dist = distance_point_circle(point1, opt2.local_disks()[id2]) + opt2.local_disks()[id2].radius;
                    min_dist = std::min(min_dist, dist);
                    // TODO : synchro matlab
                    dangling_point = false; // Discontinuity in the disk list
                }
            } // Loop over disks of obj2

            // If there is a valid dangling point, it must be on the closing point of the other floe boundary
            if (dangling_point && dangling_id == get_floe(n2).geometry().outer().size())
            {
                const segment_type segment = segment_from_id1( n2, dangling_id );
                if (segment_pos(segment, point1) <= 0)
                {
                    // Contact point-point (check sub-derivative)
                    const point_type point2 = point_from_id(n2, dangling_id);
                    const value_type dist = distance(point1, point2);
                    if (
                            dist < min_dist
                            &&  segment_pos( segment_from_id2(n1, ipt1), point2 ) >= 1
                            &&  segment_pos( segment_from_id1(n1, ipt1), point2 ) <= 0
                       )
                    {
                        min_contact = create_contact(n1, n2, point1, point2);
                        min_dist = dist;
                        // seg_id = -1; // TEST
                    }
                }
            }

            // Add contact if any
            //if (min_dist <= opt2.cdist())
            if (min_dist <= std::min( opt1.cdist(), opt2.cdist() ) )
            {
                /*
                // TEST
                if (seg_id >= 0)
                {
                    if (segment_cnt[seg_id] < 2)
                    {
                        contact_list.push_back(min_contact);
                        ++segment_cnt[seg_id];
                    }
                } else
                    contact_list.push_back(min_contact);
                */
                contact_list.push_back(min_contact);
            }

            global_min_dist = std::min( global_min_dist, min_dist );

        } // Loop over points of this disk

    } // Loop over disks of obj1

    // Add edge in graph if there is any contact
    if (contact_list.size() != 0)
    {
        add_edge(vertex(real_floe_id(n1), m_contacts), vertex(real_floe_id(n2), m_contacts), {contact_list, n1, n2}, m_contacts);
        set_dist_opt(n1, n2, std::min(opt1.cdist(), opt2.cdist()));
    }

    // Return minimal distance between the 2 floes
    return global_min_dist;
}


template <
    typename TFloe,
    typename TContact
>
void
MatlabDetector<TFloe, TContact>::clean_dist_opt()
{
    for ( auto const& edge : make_iterator_range( edges( m_contacts ) ) )
    {
        auto const& contact = m_contacts[edge];
        if (!contact.is_solved())
            set_dist_opt(contact.n1(), contact.n2(), 0);
    }
}


/*! Relative position of the projection of point on the segment
 *
 * \param segment the segment
 * \param point point to project on the segment
 */
template <
    typename TFloe,
    typename TContact
>
typename MatlabDetector<TFloe, TContact>::value_type
MatlabDetector<TFloe, TContact>::segment_pos( segment_type const& segment, point_type const& point ) const
{
    const point_type u = segment.second - segment.first;

    return floe::geometry::dot_product( u, point - segment.first ) / sum(u*u);
}

/*! Distance of a point to a segment
 *
 * \param segment the segment
 * \param point   point to project on the segment
 */
template <
    typename TFloe,
    typename TContact
>
typename MatlabDetector<TFloe, TContact>::value_type
MatlabDetector<TFloe, TContact>::segment_dist( segment_type const& segment, point_type const& point ) const
{
    const point_type u = segment.second - segment.first;
    u = { -u[1], u[0] };

    return floe::geometry::dot_product( u, point - segment.first ) / norm2(u);
}

/*! Segment from a floe beginning with specified point id
 *
 * \param n     floe id
 * \param id1   point id
 */
template <
    typename TFloe,
    typename TContact
>
typename MatlabDetector<TFloe, TContact>::segment_type
MatlabDetector<TFloe, TContact>::segment_from_id1( std::size_t n, std::size_t id1 ) const
{
    auto const& boundary = get_floe(n).geometry().outer();
    
    std::size_t id2 = id1 + 1;
    
    if (id1 >= boundary.size())    id1 -= boundary.size();
    if (id2 >= boundary.size())    id2 -= boundary.size();

    return { boundary[id1], boundary[id2] };
}

/*! Segment from a floe ending with specified point id
 *
 * \param n     floe id
 * \param id2   point id
 */
template <
    typename TFloe,
    typename TContact
>
typename MatlabDetector<TFloe, TContact>::segment_type
MatlabDetector<TFloe, TContact>::segment_from_id2( std::size_t n, std::size_t id2 ) const
{
    auto const& boundary = get_floe(n).geometry().outer();
    
    if (id2 >= boundary.size())    id2 -= boundary.size();
    std::size_t id1 = id2 + boundary.size() - 1;
    if (id1 >= boundary.size())    id1 -= boundary.size();

    return { boundary[id1], boundary[id2] };
}

/*! Point from a floe with specified id
 *
 * \param n     floe id
 * \param id    point id
 */
template <
    typename TFloe,
    typename TContact
>
typename MatlabDetector<TFloe, TContact>::point_type
MatlabDetector<TFloe, TContact>::point_from_id( std::size_t n, std::size_t id ) const
{
    auto const& boundary = get_floe(n).geometry().outer();
    
    if (id >= boundary.size())    id -= boundary.size();

    return boundary[id];
}

/*! Point for his position on a segment
 *
 * \param segment   The segment.
 * \param pos       The relative position of the point
 */
template <
    typename TFloe,
    typename TContact
>
inline
typename MatlabDetector<TFloe, TContact>::point_type
MatlabDetector<TFloe, TContact>::point_from_pos( segment_type const& segment, value_type pos ) const
{
    return (1.-pos) * segment.first + pos * segment.second;
}

}}} // namespace floe::collision::matlab

#endif // FLOE_COLLISION_MATLAB_DETECTOR_HPP
