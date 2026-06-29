/*!
 * \file generator.hpp
 * \brief Floe generator
 * \author Quentin Jouet
 */

#ifndef GENERATOR_GENERATOR_DEF_HPP
#define GENERATOR_GENERATOR_DEF_HPP

#include "floe/generator/generator.h"
#include "floe/generator/mesh_generator.hpp"

// Boost geometry
#include "floe/geometry/frame/frame_transformers.hpp"
#include <cmath>
#include <cstddef>     // std::size_t
#include <type_traits> // std::enable_if
#include <ostream>
#include "floe/arithmetic/container_operators.hpp"

// Matlab io
#include <matio.h>
// HDF5 io (new floe-shape libraries; see make_biblio_h5.py)
#include "H5Cpp.h"

#include <ctime>
#include <algorithm>
#include <random>
#include <unordered_map>
#include "floe/utils/random.hpp"

// assertion
#include <cassert>

// TEST
#include "floe/integration/gauss_legendre.hpp"
#include "floe/integration/integrate.hpp"

namespace floe { namespace generator {

using namespace types;

// CGAL typedefs
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Triangulation_vertex_base_2<K>;
using Fb = CGAL::Delaunay_mesh_face_base_2<K>;
using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, Tds>;
using Criteria = CGAL::Delaunay_mesh_size_criteria_2<CDT>;
using Vertex_handle = CDT::Vertex_handle;
using Point = CDT::Point;

template<typename T>
using scale_transformer = boost::geometry::strategy::transform::scale_transformer<T, 2, 2>;


template<typename TProblem>
void
Generator<TProblem>::generate_floe_set(std::size_t nb_floes, real_type concentration, real_type max_size, real_type min_size, std::vector<int> force_modes,
    std::vector<real_type> force_speeds)
{   
    // basic process : 
    // * floes are generated in spiral, with the biggest at the center. Convergent wind/current are applied to gather the floes  
    // * regularly, all floes within a subwindow are stopped to reduce the kinetic energy
    // * the size of this subwindow is determined by the shrinkFactor 
    // * the number of floes outside the window is looked at to decide wether it is time to exit the loop 
    // * all floes are stopped at the end 
    double shrinkFactor(0.85);
    std::cout << "Generate " << nb_floes << " floes..." << std::endl;
    std::cout << "the restitution coefficient is fixed to: " << 
    m_problem.get_lcp_manager().get_solver().get_epsilon() << std::endl;
    load_biblio_floe(m_biblio_path); // CLI --biblio (default: realistic HDF5 library)
    discretize_biblio_floe(25);
    generate_meshes();
    random_floe_group(nb_floes, max_size, min_size);
    real_type mu_static = 0;
    m_problem.get_floe_group().set_mu_static(mu_static);
    std::cout << "the ice/ice static friction coefficient is fixed to: " << mu_static << std::endl;
    assert( (m_problem.get_lcp_manager().get_solver().get_epsilon()==0 && mu_static==0) 
        && "Error: to optimize the initial configuration generation,\n"
        && "please set the restitution and the ice/ice static friction coefficients to zero!\n" );

    m_problem.create_optim_vars();
    real_type max_radius = 0;
    real_type min_radius = std::numeric_limits<real_type>::max();
    for (auto optim_ptr : m_problem.proximity_detector().data().get_optims()){
        real_type floe_global_disk_radius = optim_ptr->global_disk().radius;
        max_radius = std::max(max_radius, floe_global_disk_radius);
        min_radius = std::min(min_radius, floe_global_disk_radius);
    }
    std::cout << "Floe diameter : "
        << "max = " << max_radius * 2
        << ", min = " << min_radius * 2 << std::endl;

    real_type win_width = sqrt(m_problem.get_floe_group().total_area() / concentration);
    std::cout << "WIDTH WINDOW : " << win_width << std::endl;
    m_window = {{-win_width / 2, win_width / 2, -win_width / 2, win_width / 2}};
    m_problem.get_floe_group().set_initial_window(m_window);
    auto& physical_data = m_problem.get_dynamics_manager().get_external_forces().get_physical_data();
    // physical_data.set_window_size(win_width * 0.99, win_width * 0.99);
    physical_data.set_window_size(win_width * shrinkFactor, win_width * shrinkFactor);

    assert( (force_modes[0]==2 && force_modes[1]==0) || (force_modes[0]==0 && force_modes[1]==2) 
        && "Error : The Atmospheric or/and Oceanic currents are not suitable for initial ice pack generation.\n" 
        && "Please type:   air mode = 2 and water mode = 0 or air mode = 0 and water mode = 2.\n");
    physical_data.set_modes(force_modes[0],force_modes[1]);

    physical_data.set_speeds(force_speeds[0],force_speeds[1]);
    for (auto& floe : m_problem.get_floe_group().get_floes()){
        floe.static_floe().set_thickness(floe.static_floe().thickness() * max_size / 250);
    }
    real_type end_time = 20;
    real_type time_to_stop_floe_in_target_window = 50*end_time;//1000;
    // real_type time_to_stop_floe_in_target_window = 1000;//1000;
    bool init = true;
    
    m_problem.set_is_generator();
    do {
        m_problem.solve(end_time, 10, 10, init);
        // Trick for helping floe to stay within the target area 
        if (end_time>time_to_stop_floe_in_target_window) {
            std::cout << "It is time to stop floes within the target area." << std::endl;
            m_problem.get_floe_group().stop_floes_in_window(win_width * shrinkFactor, win_width * shrinkFactor);
            // time_to_stop_floe_in_target_window += 500; //500
            time_to_stop_floe_in_target_window += 5*end_time; //500
        }
        std::cout << " Concentration : " << m_problem.floe_concentration() << std::endl;
        end_time += 20; init = false;
        if (*m_problem.QUIT) break; // exit normally after SIGINT
    } 
    // while (m_problem.get_floe_group().kinetic_energy() != 0 && end_time < 1e6 );
    while (m_problem.get_floe_group().count_floes_outside_window(win_width*0.999, win_width*0.999) > 0 && end_time < 1e6 );
    m_problem.get_floe_group().stop_floes_in_window(win_width, win_width);
    // m_problem.get_floe_group().stop_floes_in_window(win_width, win_width); // not sure that's useful, just in case the velocities are written in the input file 

    //!< \remark    It is better do not use the while condition below since for building a big floe packs as
    //!<            an assembly of generated 2000-floe packs one need for exact inclusion within a box, otherwise
    //!<            interpenetrations may occur at the border!
    //!<            A good way for exact inclusion is to require the kinetic energy equal to zero since this is
    //!<            equivalent to floes are inside (without contact with the border) the box. 
    // } while (m_problem.get_floe_group().kinetic_energy() != 0 && end_time < 1e6 && (concentration-m_problem.floe_concentration() > 2e-3) );
    // m_problem.get_floe_group().stop_floes_in_window(win_width, win_width);
    m_problem.get_floe_group().reset_impulses();
    // set thickness to 0 for having the same initial configuration as the one generated with the input file (since thickness was not present in the input file before 2023)
    for (auto& floe : m_problem.get_floe_group().get_floes()){
        floe.static_floe().set_thickness(0);
    }
}

template<typename TProblem>
void
Generator<TProblem>::random_floe_group(std::size_t n, real_type max_size, real_type min_size)
{
    auto& list_floes = m_problem.get_floe_group().get_floes();
    // Floe-size distribution selected by --sizerep (see set_size_rep)
    std::vector<real_type> sizes;
    switch (m_size_rep)
    {
        case 2:  sizes = two_sizes_repartition(n, max_size); break;
        case 3:  sizes = random_size_repartition(n, max_size); break;
        default: sizes = exp_size_repartition(n, max_size, min_size); break;
    }
    std::vector<double>::iterator result = std::min_element(std::begin(sizes), std::end(sizes));
    std::cout << "min diameter: " << *result << "\n";
    auto min_s = *result;
    result = std::max_element(std::begin(sizes), std::end(sizes));
    std::cout << "max diameter: " << *result << "\n";
    std::cout << "The size ratio is: " << *result/min_s << std::endl;
    // Largest circumscribed radius (about the centroid) over the library, per unit nominal size: shapes in
    // m_biblio_floe_h are already recentred on their mass centre (generate_meshes) and normalized ~unit, so
    // a floe scaled by sizes[i] has true circumradius <= circum_factor * sizes[i]. scattered_distribution
    // bounds every floe by this max (each draws a random shape only afterwards) to stay overlap-free.
    real_type circum_factor = 1;
    for (auto const& s : m_biblio_floe_h)
        for (auto const& p : s.outer())
            circum_factor = std::max(circum_factor, std::sqrt(p.x * p.x + p.y * p.y));
    // Initial placement layout (CLI --distrib): 0 = uniform scatter (default), 1 = legacy spiral.
    auto centers = (m_distrib == 1) ? spiral_distribution(sizes, max_size)
                                    : scattered_distribution(sizes, max_size, circum_factor);
    auto generator = floe::random::get_uniquely_seeded_generator();
    std::uniform_int_distribution<int> distribution(0, m_biblio_size - 1);
    // unused pick to avoid first choice often being the same
    // (for quasi-simultaneous execs, even with different seeds)
    distribution(generator);
    list_floes.resize(n);

    for (std::size_t i = 0; i < n; i++)
    {
        int idx = distribution(generator);
        // int idx = 0;
        auto& base_shape = m_biblio_floe_h[idx];
        auto mesh = m_biblio_floe_h_meshes[idx];
        polygon_type shape;
        auto size = sizes[i];
        geometry::transform( base_shape, shape, scale_transformer<real_type>{ size } );
        geometry::transform( mesh, mesh, scale_transformer<real_type>{ size } );

        // Create Kinematic floe
        auto& floe = list_floes[i];
        // link static floe
        floe.attach_static_floe_ptr(std::unique_ptr<static_floe_type>(new static_floe_type()));
        
        auto& static_floe = floe.static_floe();

        // Import boundary
        std::unique_ptr<typename floe_type::geometry_type> geometry(new typename floe_type::geometry_type(shape));
        static_floe.attach_geometry_ptr(std::move(geometry));
        
        // Import mesh
        floe.static_floe().set_mesh(mesh);
        // Done.

        // Set space-time state
        typename floe_type::state_type state {{0, 0}, 0, {0, 0}, 0, {0, 0}};
        state.pos = centers[i];
        floe.set_state( state );
    }

}


template<typename TProblem>
std::vector<real_type>
Generator<TProblem>::random_size_repartition(std::size_t n, real_type R_max)
{
    std::vector<real_type> v;
    std::default_random_engine generator;
    std::exponential_distribution<double> distribution(5);
    real_type R_min = 0;
    for (std::size_t i = 0; i < n; i++)
    {
        real_type R = R_min + std::min(distribution(generator), 1.) * (R_max - R_min);
        v.push_back(R);
    }
    return v;
}

template<typename TProblem>
std::vector<real_type>
Generator<TProblem>::exp_size_repartition(std::size_t n, real_type R_max, real_type R_min)
{
    std::random_device rd;
    std::mt19937 g(rd());

    std::cout << "The exponent of the power law is: " << m_alpha << " and the floe number per size is: " << m_nbfpersize << std::endl;
    std::vector<real_type> v;
    // allow to generate more than one floe per size categories! 
    // Useful for making easier the init. config. generation!
    // Yet, warning, since the exponential distribution of size is no longer corresponding to m_alpha!! 
    for (std::size_t i = 1; i <= n/m_nbfpersize; i++)
    {
        real_type R =  std::max(R_max * exp( (- 1 /  m_alpha) * log(i) ), R_min);
        for (int j=0; j<m_nbfpersize; ++j){
            v.push_back(R);
        }
    }
    // std::srand ( unsigned ( std::time(0) ) ); // seed for not having pseudorandom
    std::shuffle ( v.begin() + 1, v.end(), g); // first floe(biggest) stays first (-> initially at center)
    // std::random_shuffle ( v.begin() + 1, v.end() ); // first floe(biggest) stays first (-> initially at center)
    return v;
}

/*! Two-size distribution: each floe is R_max or R_max/ratio (~half each). Useful for binary-size packs.
 *  The first floe is forced to R_max so the biggest stays at the spiral center, like exp_size_repartition. */
template<typename TProblem>
std::vector<real_type>
Generator<TProblem>::two_sizes_repartition(std::size_t n, real_type R_max, real_type ratio)
{
    std::random_device rd;
    std::mt19937 g(rd());
    std::cout << "Two-size distribution: " << R_max << " and " << R_max / ratio << " (~half each)" << std::endl;
    std::vector<real_type> v;
    std::bernoulli_distribution coin(0.5);
    for (std::size_t i = 0; i < n; ++i)
        v.push_back(coin(g) ? R_max : R_max / ratio);
    if (!v.empty()) v[0] = R_max;                 // biggest first -> initially at center
    std::shuffle(v.begin() + 1, v.end(), g);
    return v;
}

/*! Uniform random scatter of overlap-free floes (stage A of pack prep). Unlike the spiral, it has no
 *  radial bias, so after the convergent compaction there is no central hole nor edge ring to fix.
 *  Algorithm: random sequential placement, biggest floes first (hardest to fit), rejecting any position
 *  whose bounding disk overlaps an already-placed one; the search region (sized for a comfortable SEED
 *  concentration) enlarges when a floe can't be placed, so late/over-crowded floes simply spill outside
 *  the window — the convergent forcing (mode 2) then gathers them in. A spatial hash with per-floe
 *  register/query radii keeps it O(N) despite the 100-200x size dispersity. */
template<typename TProblem>
std::vector<point_type>
Generator<TProblem>::scattered_distribution(std::vector<real_type> const& sizes, real_type Rmax,
    real_type circum_factor)
{
    const std::size_t n = sizes.size();
    std::vector<point_type> centers(n);
    if (n == 0) return centers;

    // Bound each floe by its CIRCUMSCRIBED disk (circum_factor * nominal size): the nominal size is the
    // 1/Rmin normalization, but the real recentred shape extends up to circum_factor times that, so a disk
    // of radius = size alone underestimates the floe and lets neighbours interpenetrate (seen with
    // --nbfpersize > 1, where bigger floes make the overlap large enough to fire an INTER at t=0).
    real_type total_disk_area = 0;
    for (auto s : sizes) total_disk_area += M_PI * (circum_factor * s) * (circum_factor * s);
    const real_type SEED_CONC = 0.45; // bounding-disk fill of the seed region (real floe fill is lower)
    real_type half = 0.5 * std::sqrt(total_disk_area / SEED_CONC); // half-side of the square seed region

    const real_type cell = std::max(circum_factor * Rmax / 4, std::numeric_limits<real_type>::min());
    const real_type GAP = 1.02; // small no-overlap safety margin
    std::unordered_map<long long, std::vector<std::size_t>> grid;
    auto key = [](long long ix, long long iy) { return (ix + 2000000LL) * 4000000LL + (iy + 2000000LL); };
    std::vector<point_type> pos(n);
    std::vector<real_type> rad(n);

    auto cells_radius = [&](real_type r) { return (long long)std::ceil(r / cell) + 1; };
    auto overlaps = [&](point_type p, real_type r) -> bool {
        const long long c = cells_radius(r);
        const long long cx = (long long)std::floor(p.x / cell), cy = (long long)std::floor(p.y / cell);
        for (long long ix = cx - c; ix <= cx + c; ++ix)
            for (long long iy = cy - c; iy <= cy + c; ++iy) {
                auto it = grid.find(key(ix, iy));
                if (it == grid.end()) continue;
                for (auto j : it->second) {
                    const real_type dx = p.x - pos[j].x, dy = p.y - pos[j].y;
                    const real_type rr = (r + rad[j]) * GAP;
                    if (dx * dx + dy * dy < rr * rr) return true;
                }
            }
        return false;
    };
    auto place = [&](std::size_t idx, point_type p, real_type r) {
        pos[idx] = p; rad[idx] = r; centers[idx] = p;
        const long long c = cells_radius(r);
        const long long cx = (long long)std::floor(p.x / cell), cy = (long long)std::floor(p.y / cell);
        for (long long ix = cx - c; ix <= cx + c; ++ix)
            for (long long iy = cy - c; iy <= cy + c; ++iy)
                grid[key(ix, iy)].push_back(idx);
    };

    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return sizes[a] > sizes[b]; });

    auto gen = floe::random::get_uniquely_seeded_generator();
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t i = order[k];
        const real_type r = circum_factor * sizes[i]; // circumscribed bounding radius (see above)
        real_type reach = half;
        for (int attempt = 0; ; ++attempt) {
            std::uniform_real_distribution<real_type> u(-reach, reach);
            point_type p{ u(gen), u(gen) };
            if (!overlaps(p, r)) { place(i, p, r); break; }
            if ((attempt + 1) % 40 == 0) reach *= 1.08; // struggling: enlarge (floe spills outside window)
        }
    }
    return centers;
}

template<typename TProblem>
std::vector<point_type>
Generator<TProblem>::spiral_distribution(std::vector<real_type> const& size_distribution, real_type Rmax)
{
    std::vector<point_type> v;
    real_type R = 0, theta = 0, R_max = *std::max_element(size_distribution.begin(), size_distribution.end());

    std::default_random_engine random_gen;
    std::uniform_real_distribution<real_type> uniform_distrib;
    // for (real_type r : size_distribution)
    for (auto it = size_distribution.begin(); it != size_distribution.end(); ++it)
    {
        real_type d_theta, d_R, r = *it;
        if (R < R_max)
            d_theta = r * 1.2 / R_max, d_R = r * 2;
        else
            d_theta = r * 1.2 / R, d_R = (Rmax / M_PI + 1) * 1.5 * (r / R);

        if (it == size_distribution.begin())
        {
            v.push_back(point_type{R * cos(theta), R * sin(theta)});
            R += d_R;
        } else {
            theta += d_theta; R += d_R;
            real_type Rf = R + uniform_distrib(random_gen) * (R_max - r) * 0.8;
            v.push_back(point_type{Rf * cos(theta), Rf * sin(theta)});
            theta += d_theta; R += d_R;
        }
        
    }
    return v;
}

template<typename TProblem>
void Generator<TProblem>::load_biblio_floe(std::string filename)
{
    // Dispatch on extension: ".h5" -> new HDF5 library, otherwise legacy matio (.mat).
    if (filename.size() >= 3 && filename.compare(filename.size() - 3, 3, ".h5") == 0)
        load_biblio_floe_h5(filename);
    else
        load_biblio_floe_mat(filename);
}

/*! HDF5 floe-shape library loader (schema written by pack_creator/make_biblio_h5.py):
 *  group "shapes" with one (N,2) dataset per shape ("0".."M-1"), a "Rmin" (M) vector and a "Cmin"
 *  (M,2). Each shape is normalized to radius 1 (scaled by 1/Rmin), exactly like the matio path. */
template<typename TProblem>
void Generator<TProblem>::load_biblio_floe_h5(std::string filename)
{
    using namespace H5;
    H5File file(filename, H5F_ACC_RDONLY);

    // radii (M) -> also gives the number of shapes
    DataSet rmin_ds = file.openDataSet("Rmin");
    hsize_t rdim[1] = {0};
    rmin_ds.getSpace().getSimpleExtentDims(rdim, NULL);
    m_biblio_size = rdim[0];
    std::vector<real_type> rmin(m_biblio_size);
    rmin_ds.read(rmin.data(), PredType::NATIVE_DOUBLE);

    Group shapes = file.openGroup("shapes");
    for (std::size_t i = 0; i < m_biblio_size; ++i)
    {
        DataSet ds = shapes.openDataSet(std::to_string(i));
        hsize_t dims[2] = {0, 0};
        ds.getSpace().getSimpleExtentDims(dims, NULL); // dims[0] = N points, dims[1] = 2
        const std::size_t n = dims[0];
        std::vector<real_type> buf(n * 2);
        ds.read(buf.data(), PredType::NATIVE_DOUBLE);  // row-major: [x0,y0,x1,y1,...]

        multi_point_type shape;
        for (std::size_t j = 0; j < n; ++j)
            shape.push_back(point_type{buf[2 * j], buf[2 * j + 1]});

        const real_type radius = rmin[i];
        geometry::transform(shape, shape, scale_transformer<real_type>{1 / radius});
        m_biblio_floe.push_back(shape);
    }
}

template<typename TProblem>
void Generator<TProblem>::load_biblio_floe_mat(std::string filename)
{
    mat_t *matfp;

    // Opening file
    matfp = Mat_Open( filename.c_str(), MAT_ACC_RDONLY );
    if ( matfp == nullptr )
    {
        throw std::ios_base::failure("Error opening MAT file \"" + filename + "\"");
    }

    matvar_t *shape_list, *Rmin, *Cmin;

    shape_list = Mat_VarRead(matfp,"G");
    Rmin = Mat_VarRead(matfp,"Rmin"); // surrounding disks radius
    Cmin = Mat_VarRead(matfp,"Cmin"); // surrounding disks centers

    m_biblio_size = shape_list->dims[0];

    // checking vars
    for (auto* matvar : {shape_list, Rmin, Cmin})
    {
        if ( NULL == matvar ) {
            fprintf(stderr,"Variable not found, or error reading MAT file\n");
        }
    }

    // Get the cells
    matvar_t **cells = (matvar_t **)shape_list->data;

    // Import shapes
    for (std::size_t i = 0; i < m_biblio_size; i++)
    {   
        // Read geometry
        auto cell = cells[i];
        multi_point_type shape;
        for (std::size_t j = 0; j < cell->dims[0]; j++)
        {
            shape.push_back(point_type{static_cast<real_type*>(cell->data)[j], static_cast<real_type*>(cell->data)[cell->dims[0] + j]});
        }

        // Center and normalize geometry
        real_type radius{static_cast<real_type*>(Rmin->data)[i]};
        point_type center{static_cast<real_type*>(Cmin->data)[i], static_cast<real_type*>(Cmin->data)[m_biblio_size + i]};
        // geometry::transform( shape, shape, geometry::frame::transformer( typename floe_type::frame_type{-center, 0} )); // done when creating mesh
        geometry::transform( shape, shape, scale_transformer<real_type>{ 1 / radius } );

        // Save geometry
        m_biblio_floe.push_back(shape);
    }

    // freeing memory
    Mat_VarFree(shape_list);

    Mat_Close(matfp);
}

template<typename TProblem>
void Generator<TProblem>::discretize_biblio_floe(std::size_t n)
{
    for (auto shape : m_biblio_floe)
    {
        // TODO find a lib to do that properly (equivalent DecimatePoly.m)
        polygon_type shape_h;
        std::size_t nb_pts = std::min(shape.size(),  n);
        for (std::size_t i = 0; i < nb_pts; i++)
            shape_h.outer().push_back(shape[(shape.size() * i) / nb_pts]);
        // Save the simplified shape
        m_biblio_floe_h.push_back(shape_h);
    }
}


template<typename TProblem>
void Generator<TProblem>::generate_meshes()
{
    for (auto& shape : m_biblio_floe_h)
    {
        auto mesh = generate_mesh_for_shape<polygon_type, mesh_type>(shape);
        // Center mesh and shape on floe's center of mass
        using integration_strategy = floe::integration::RefGaussLegendre<real_type,2,2>;
        auto mass_center = floe::integration::integrate(
            [] (real_type x, real_type y) { return point_type{x, y}; },
            mesh, integration_strategy()
        ) / floe::integration::integrate(
            [] (real_type x, real_type y) { return 1.; },
            mesh,integration_strategy()
        );
        // std::cout << "MC : " << mass_center << std::endl;
        polygon_type shape_cpy = shape;
        geometry::transform( shape_cpy, shape, geometry::frame::transformer( typename floe_type::frame_type{-mass_center, 0} ));
        mesh_type mesh_cpy = mesh;
        geometry::transform( mesh_cpy, mesh, geometry::frame::transformer( typename floe_type::frame_type{-mass_center, 0} ));
        // Save mesh
        m_biblio_floe_h_meshes.push_back(mesh);
    }
}


}} // namespace floe::generator


#endif // GENERATOR_GENERATOR_DEF_HPP
