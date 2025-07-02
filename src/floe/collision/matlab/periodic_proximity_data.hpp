/*!
 * \file floe/collision/matlab/periodic_proximity_data.hpp
 * \brief 
 * \author Quentin Jouet
 */

#ifndef FLOE_COLLISION_PERIODIC_PROXIMITY_DATA
#define FLOE_COLLISION_PERIODIC_PROXIMITY_DATA


#include "floe/collision/matlab/proximity_data.hpp"
#include "floe/collision/matlab/ghost_optimized_floe.hpp"
#include "floe/floes/ghost_floe.hpp"


namespace floe { namespace collision { namespace matlab
{

namespace ublas = boost::numeric::ublas;


template <
    typename TFloeGroup,
    typename TOptim
>
class PeriodicProximityData : public ProximityData<TFloeGroup, TOptim>
{

public:

    using base_class = ProximityData<TFloeGroup, TOptim>;
    using real_type = typename base_class::real_type;
    using floe_type = typename base_class::floe_type;
    using point_type = typename floe_type::point_type;
    using optim_type = TOptim;
    using floe_interface_type = typename floe_type::floe_interface_type;
    using optim_interface_type = typename optim_type::optim_interface_type;

    using ghost_floe_type = floes::GhostFloe<floe_type>;
    using ghost_optim_type = GhostOptimizedFloe<optim_type>;

    //! Default constructor
    PeriodicProximityData() : base_class() { }

    //!Empty floe and optim lists
    virtual void reset() override {
        base_class::reset();
        m_ghost_floes.clear();
        m_ghost_optims.clear();
    }

    void add_ghost(std::size_t floe_id, point_type translation){
        m_ghost_floes.push_back(ghost_floe_type{ this->get_floe(floe_id), translation, floe_id});
        m_ghost_optims.push_back(ghost_optim_type{ this->get_optim(floe_id), translation });
    }

    inline std::size_t nb_ghosts() const { return m_ghost_floes.size(); }

    inline virtual void set_dist_secu(std::size_t n1, std::size_t n2, real_type val) override {
        (n2 >= this->nb_floes()) ? this->m_dist_secu(n1, n2) = val : this->m_dist_secu(n1, n2) = this->m_dist_secu(n2, n1) = val;
    }
    inline virtual void set_indic(std::size_t n1, std::size_t n2, short val) override {
        (n2 >= this->nb_floes()) ? this->m_indic(n1, n2) = val : this->m_indic(n1, n2) = this->m_indic(n2, n1) = val;
    }
    virtual void set_dist_opt(std::size_t n1, std::size_t n2, real_type val) override;

    virtual floe_interface_type const& get_floe_itf(std::size_t n) const override;
    virtual optim_interface_type const& get_optim_itf(std::size_t n) const override;
    inline ghost_floe_type const& get_ghost_floe(std::size_t n) const { return m_ghost_floes[n]; }

    //! returns real floe id if param id refers to a ghost
    virtual std::size_t real_floe_id(std::size_t n) const override{
        std::size_t N { this->nb_floes() };
        if (n >= N)
            return m_ghost_floes[n - N].m_original_id;
        else
            return n;
    }

private:
    std::vector<ghost_floe_type>     m_ghost_floes; //!< Ghost floes list.
    std::vector<ghost_optim_type>    m_ghost_optims; //!< Ghost Optimized floes list.

};


template <
    typename TFloeGroup,
    typename TOptim
>
typename PeriodicProximityData<TFloeGroup, TOptim>::floe_interface_type const&
PeriodicProximityData<TFloeGroup, TOptim>
::get_floe_itf(std::size_t n) const
{
    const std::size_t N { this->nb_floes() };
    if (n < N)
        return this->get_floes()[n];
    else
        return m_ghost_floes[n - N];
}


template <
    typename TFloeGroup,
    typename TOptim
>
typename PeriodicProximityData<TFloeGroup, TOptim>::optim_interface_type const&
PeriodicProximityData<TFloeGroup, TOptim>
::get_optim_itf(std::size_t n) const
{
    const std::size_t N { this->nb_floes() };
    if (n < N)
        return base_class::get_optim_itf(n);
    else
        return m_ghost_optims[n - N];
}

template <
    typename TFloeGroup,
    typename TOptim
>
void
PeriodicProximityData<TFloeGroup, TOptim>
::set_dist_opt(
    std::size_t n1, std::size_t n2, real_type val)
{
    const std::size_t N { this->nb_floes() };
    if (n2 >= N){
        this->m_dist_opt(n1, n2) = val;
    } else if (n1 >= N) {
        this->m_dist_opt(n2, n1) = val;
    } else {
        this->m_dist_opt(n1, n2) = this->m_dist_opt(n2, n1) = val;
    }
}


}}} // namespace floe::collision::matlab

#endif // FLOE_COLLISION_PERIODIC_PROXIMITY_DATA
