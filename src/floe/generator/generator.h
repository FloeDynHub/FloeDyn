/*!
 * \file generator/generator.h
 * \brief Floe generator
 * \author Quentin Jouet
 */

#ifndef GENERATOR_GENERATOR_HPP
#define GENERATOR_GENERATOR_HPP

// Boost geometry
#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/point.hpp"   
#include "floe/geometry/geometries/multi_point.hpp"

// Floes
#include "floe/floes/kinematic_floe.hpp"

#include <vector>
#include <atomic>
#include <iostream>


namespace floe { namespace generator {

template<typename TProblem>
class Generator
{
public:
    using real_type = typename TProblem::real_type;
    using floe_type = typename TProblem::floe_group_type::floe_type;
    using point_type = typename floe_type::point_type;
    using multi_point_type = floe::geometry::MultiPoint<point_type>;
    using polygon_type = typename floe_type::geometry_type;
    using static_floe_type = typename floe_type::static_floe_type;
    using mesh_type = typename floe_type::mesh_type;

    Generator(real_type alpha,int nbfpersize) : m_problem{0.0, 0}, m_alpha{alpha}, m_nbfpersize{nbfpersize} {}

    //! Generate floe set with given number of floe and concentration
    void generate_floe_set(std::size_t number, real_type concentration, real_type max_size, 
        std::vector<int> force_modes, std::vector<real_type> force_speeds);
    typename TProblem::floe_group_type& get_floe_group() { return m_problem.get_floe_group(); }
    // std::array<real_type, 4> get_window() const { return m_window; }
    // real_type window_area() const { return (m_window[1] - m_window[0]) * (m_window[3] - m_window[2]); }
    void set_exit_signal(std::atomic<bool>* QUIT){ m_problem.QUIT = QUIT; }

    inline void set_frac_dim(real_type alpha) { m_alpha = alpha; }
    inline void set_nb_floe_per_size(int nb) { m_nbfpersize = nb; }

private:
    TProblem m_problem;
    std::array<real_type, 4> m_window;
    //! Floe size repartition (exponential)
    std::vector<real_type> random_size_repartition(std::size_t n, real_type R_max);
    std::vector<real_type> exp_size_repartition(std::size_t n, real_type R_max);
    //! Random floe group
    void random_floe_group(std::size_t n, real_type max_size);
    //! Spiral dispatcher
    std::vector<point_type> spiral_distribution(std::vector<real_type> const& size_distribution, real_type Rmax);

    void load_biblio_floe(std::string filename);
    void discretize_biblio_floe(std::size_t n);
    void generate_meshes();

    std::size_t m_biblio_size;
    std::vector<multi_point_type> m_biblio_floe;
    std::vector<polygon_type> m_biblio_floe_h;
    std::vector<mesh_type> m_biblio_floe_h_meshes;
    std::vector<static_floe_type> m_static_floe_list;

    real_type   m_alpha; //!< fractal dimension \f$ 1.2 <= \alpha <= 3 \f$ 
    int         m_nbfpersize;
};


}} // namespace floe::generator


#endif // GENERATOR_GENERATOR_HPP
