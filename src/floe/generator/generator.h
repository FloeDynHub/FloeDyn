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


namespace floe { namespace generator {

template<typename TProblem>
class Generator
{
public:
    using value_type = typename TProblem::value_type;
    using floe_type = typename TProblem::floe_group_type::floe_type;
    using point_type = typename floe_type::point_type;
    using multi_point_type = floe::geometry::MultiPoint<point_type>;
    using polygon_type = typename floe_type::geometry_type;
    using static_floe_type = typename floe_type::static_floe_type;
    using mesh_type = typename floe_type::mesh_type;

    Generator() : m_problem{} {}

    //! Generate floe set with given number of floe and concentration
    void generate_floe_set(std::size_t number, value_type concentration);
    typename TProblem::floe_group_type& get_floe_group() { return m_problem.get_floe_group(); }
    std::array<value_type, 4> get_window() const { return m_window; }
    value_type window_area() const { return (m_window[1] - m_window[0]) * (m_window[3] - m_window[2]); }
    void set_exit_signal(std::atomic<bool>* QUIT){ m_problem.QUIT = QUIT; }

private:
    TProblem m_problem;
    std::array<value_type, 4> m_window;
    //! Floe size repartition (exponential)
    std::vector<value_type> random_size_repartition(std::size_t n, value_type R_max);
    std::vector<value_type> exp_size_repartition(std::size_t n, value_type R_max);
    //! Random floe group
    void random_floe_group(std::size_t n);
    //! Spiral dispatcher
    std::vector<point_type> spiral_distribution(std::vector<value_type> const& size_distribution, value_type Rmax);

    void load_biblio_floe(std::string filename);
    void discretize_biblio_floe(std::size_t n);
    void generate_meshes();

    std::size_t m_biblio_size;
    std::vector<multi_point_type> m_biblio_floe;
    std::vector<polygon_type> m_biblio_floe_h;
    std::vector<mesh_type> m_biblio_floe_h_meshes;
    std::vector<static_floe_type> m_static_floe_list;
};


}} // namespace floe::generator


#endif // GENERATOR_GENERATOR_HPP
