/*!
 * \file floe/collision/matlab/proximity_data.hpp
 * \brief 
 * \author Quentin Jouet
 */

#ifndef FLOE_COLLISION_PROXYMITY_DATA
#define FLOE_COLLISION_PROXYMITY_DATA


#include <iostream> // DEBUG
#include <vector>

// uBlas
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>


namespace floe { namespace collision { namespace matlab
{

namespace ublas = boost::numeric::ublas;


template <
    typename TFloeGroup,
    typename TOptim
>
class ProximityData
{

public:

    using floe_group_type = TFloeGroup;
    using floe_type = typename floe_group_type::floe_type;
    using real_type = typename floe_type::real_type;
    using optim_type = TOptim;
    using dist_matrix_type = ublas::matrix<real_type>; //!< Type of distance matrix
    using indic_matrix_type = ublas::matrix<short>; //!< Type of indicator matrix
    using floe_interface_type = typename floe_type::floe_interface_type;
    using optim_interface_type = typename optim_type::optim_interface_type;

    //! Default constructor
    ProximityData() : /*m_floes{},*/ m_floe_group{nullptr}, m_optims{}, m_indic{0,0}, m_dist_secu{0,0}, m_dist_opt{0,0}, m_interpenetration{false} {}

    //!Empty floe and optim lists
    virtual void reset() { /*m_floes.clear();*/ m_optims.clear(); }

    inline void resize(std::size_t N1, std::size_t N2) { 
        m_indic.resize(N1, N2);
        m_dist_opt = ublas::scalar_matrix<real_type>(N1, N2, 0);
        m_dist_secu = ublas::scalar_matrix<real_type>(N1, N2, 0);
    }

    // virtual void push_back( floe_type * floe_ptr )
    // {
    //     m_floes.push_back(floe_ptr);
    //     m_optims.push_back( new optim_type{*floe_ptr} );
    // }

    virtual void set_floe_group(floe_group_type const& floe_group){
        m_floe_group = &floe_group;
        m_optims.clear();
        for (auto const& floe : get_floes())
            m_optims.push_back( new optim_type{floe} );
    }

    virtual void update_optim_vars(){
        // Add missing optims (new floes in m_floe_group)
        for (std::size_t i = m_optims.size() ; i < m_floe_group->get_floes().absolute_size() ; ++i)
            m_optims.push_back( new optim_type{m_floe_group->get_floes()(i)} ); // get_floes()(i) is the filter_off getter
    }

    inline std::size_t size1() const { return m_indic.size1(); }
    inline std::size_t size2() const { return m_indic.size2(); }
    inline std::size_t nb_floes() const { return get_floes().size(); }

    inline real_type get_dist_secu(std::size_t n1, std::size_t n2) const { return m_dist_secu(n1, n2); }
    inline short get_indic(std::size_t n1, std::size_t n2) const { return m_indic(n1, n2); }
    inline real_type get_dist_opt(std::size_t n1, std::size_t n2) const { return m_dist_opt(n1, n2); }

    inline virtual void set_dist_secu(std::size_t n1, std::size_t n2, real_type val) { m_dist_secu(n1, n2) = m_dist_secu(n2, n1) = val; }
    inline virtual void set_indic(std::size_t n1, std::size_t n2, short val) { m_indic(n1, n2) = m_indic(n2, n1) = val; }
    inline virtual void set_dist_opt(std::size_t n1, std::size_t n2, real_type val) { m_dist_opt(n1, n2) = m_dist_opt(n2, n1) = val; }

    //! Container accessors
    // inline std::vector<floe_type const*> const& get_floes() const { return m_floes; }
    inline typename floe_group_type::floe_list_type const& get_floes() const { return m_floe_group->get_floes(); }
    inline std::vector<optim_type *> const& get_optims() const { return m_optims; }
    inline virtual floe_interface_type const& get_floe_itf(std::size_t n) const { return get_floe(n); }
    inline virtual optim_interface_type const& get_optim_itf(std::size_t n) const { return get_optim(n); }
    // inline floe_type const& get_floe(std::size_t n) const { return *(m_floes[n]); }
    inline floe_type const& get_floe(std::size_t n) const { return get_floes()[n]; }
    // inline optim_type& get_optim(std::size_t n) const { return *(m_optims[n]); }
    inline optim_type& get_optim(std::size_t n) const { return *(m_optims[m_floe_group->absolute_id(n)]); }
    inline virtual std::size_t real_floe_id(std::size_t n) const { return n; }

    //! interpenetration bool accessor (true if no interpenetration)
    inline bool interpenetration() const { return m_interpenetration; }
    inline void interpenetration(bool b) { m_interpenetration = b; }

protected:

    floe_group_type const* m_floe_group;
    // std::vector<floe_type const*>   m_floes; //!< Floes list.
    std::vector<optim_type*>   m_optims; //!< Optimization datas list.
    indic_matrix_type m_indic; //!< Indicator of collision (0=far away, 1=close, 2=contact)
    dist_matrix_type m_dist_secu; //!< Security distance
    dist_matrix_type m_dist_opt; //!< Optimial distance
    bool m_interpenetration; //! Floe interpenetration

};


}}} // namespace floe::collision::matlab

#endif // FLOE_COLLISION_PROXYMITY_DATA
