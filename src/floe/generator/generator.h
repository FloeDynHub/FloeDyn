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
#include "../product/config/config.hpp"

#include <vector>
#include <atomic>
#include <iostream>



namespace floe { namespace generator {

using namespace types;

template<typename TProblem>
class Generator
{
public:
    using floe_type = types::floe_type;
    using multi_point_type = floe::geometry::MultiPoint<point_type>;
    using polygon_type = typename floe_type::geometry_type;
    using static_floe_type = typename floe_type::static_floe_type;
    using mesh_type = typename floe_type::mesh_type;

    Generator(real_type alpha,int nbfpersize) : m_problem{0.0, 0}, m_alpha{alpha}, m_nbfpersize{nbfpersize} {}

    //! Generate floe set with given number of floe and concentration
    void generate_floe_set(std::size_t number, real_type concentration, real_type max_size, real_type min_size,
        std::vector<int> force_modes, std::vector<real_type> force_speeds);
    floe_group_type& get_floe_group() { return m_problem.get_floe_group(); }
    //! Access the generator's own collision manager (distinct from the outer Problem's): lets the
    //! caller enable OPTIMJAM on the generation loop itself (see simu_runner — needs --jam_unanchored
    //! since the generated pack has no obstacles).
    auto& get_lcp_manager() { return m_problem.get_lcp_manager(); }
    //! Access the generator's own out_manager: lets the caller name the GENERATION output file (e.g.
    //! "<output>_gen" so a UI can follow the generation run — see simu_runner --output handling).
    auto& get_out_manager() { return m_problem.get_out_manager(); }
    // std::array<real_type, 4> get_window() const { return m_window; }
    // real_type window_area() const { return (m_window[1] - m_window[0]) * (m_window[3] - m_window[2]); }
    void set_exit_signal(std::atomic<bool>* QUIT){ m_problem.QUIT = QUIT; }

    inline void set_frac_dim(real_type alpha) { m_alpha = alpha; }
    inline void set_nb_floe_per_size(int nb) { m_nbfpersize = nb; }
    //! Floe-shape library used to draw each generated floe's shape (CLI --biblio). A path ending in
    //! ".h5" is read by the HDF5 loader, otherwise by the legacy matio (.mat) loader. Default is the
    //! realistic HDF5 library; pass io/library/biblio_circle.h5 for circular floes.
    inline void set_biblio_path(std::string const& path) { if (!path.empty()) m_biblio_path = path; }
    //! Floe-size distribution used by the generator (CLI --sizerep): 1 = exp_size_repartition (power law,
    //! default), 2 = two_sizes_repartition (R_max and R_max/1.4, ~half each), 3 = random_size_repartition.
    inline void set_size_rep(int s) { if (s > 0) m_size_rep = s; }
    //! Initial floe-placement layout (CLI --distrib): 0 = scattered_distribution (uniform random, default),
    //! 1 = spiral_distribution (legacy, biggest at centre). See random_floe_group dispatch.
    inline void set_distrib(int d) { if (d >= 0) m_distrib = d; }

private:
    TProblem m_problem;
    std::array<real_type, 4> m_window;
    //! Floe-size distributions (selected by m_size_rep / --sizerep)
    std::vector<real_type> random_size_repartition(std::size_t n, real_type R_max);
    std::vector<real_type> exp_size_repartition(std::size_t n, real_type R_max, real_type R_min);
    std::vector<real_type> two_sizes_repartition(std::size_t n, real_type R_max, real_type ratio = 1.4);
    //! Random floe group
    void random_floe_group(std::size_t n, real_type max_size, real_type min_size);
    //! Spiral dispatcher (biggest at centre, spiralling out)
    std::vector<point_type> spiral_distribution(std::vector<real_type> const& size_distribution, real_type Rmax);
    //! Uniform random (overlap-free) scatter: no spiral => no central hole / ring after compaction.
    //! circum_factor = library's max circumscribed-radius / nominal-size, so the bounding disks cover the
    //! real (off-round) floe extent whatever shape each floe later draws.
    std::vector<point_type> scattered_distribution(std::vector<real_type> const& size_distribution, real_type Rmax,
        real_type circum_factor = 1);

    void load_biblio_floe(std::string filename);     //!< dispatches on extension (.h5 vs .mat)
    void load_biblio_floe_h5(std::string filename);  //!< HDF5 library loader (see make_biblio_h5.py)
    void load_biblio_floe_mat(std::string filename); //!< legacy matio (.mat) library loader
    void discretize_biblio_floe(std::size_t n);
    void generate_meshes();

    std::size_t m_biblio_size;
    std::vector<multi_point_type> m_biblio_floe;
    std::vector<polygon_type> m_biblio_floe_h;
    std::vector<mesh_type> m_biblio_floe_h_meshes;
    std::vector<static_floe_type> m_static_floe_list;

    real_type   m_alpha; //!< fractal dimension \f$ 1.2 <= \alpha <= 3 \f$
    int         m_nbfpersize;
    std::string m_biblio_path{"io/library/biblio_realistic.h5"}; //!< floe-shape library (CLI --biblio)
    int         m_size_rep{1}; //!< floe-size distribution selector (CLI --sizerep; see set_size_rep)
    int         m_distrib{0};  //!< initial placement layout (CLI --distrib; 0=scattered, 1=spiral)
};


}} // namespace floe::generator


#endif // GENERATOR_GENERATOR_HPP
