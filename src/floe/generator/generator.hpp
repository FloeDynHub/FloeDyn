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

#include <ctime>
#include <algorithm>
#include <random>
#include "floe/utils/random.hpp"


// TEST
#include "floe/integration/gauss_legendre.hpp"
#include "floe/integration/integrate.hpp"

namespace floe { namespace generator {

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
Generator<TProblem>::generate_floe_set(std::size_t nb_floes, real_type concentration, real_type max_size)
{   
    std::cout << "Generate " << nb_floes << " floes..." << std::endl;
    load_biblio_floe("io/inputs/Biblio_Floes.mat");
    discretize_biblio_floe(25);
    generate_meshes();
    random_floe_group(nb_floes, max_size);
    m_problem.get_floe_group().set_mu_static(0);

    real_type win_width = sqrt(m_problem.get_floe_group().total_area() / concentration);
    std::cout << "WIDTH WINDOW : " << win_width << std::endl;
    m_window = {{-win_width / 2, win_width / 2, -win_width / 2, win_width / 2}};
    m_problem.get_floe_group().set_initial_window(m_window);
    auto& physical_data = m_problem.get_dynamics_manager().get_external_forces().get_physical_data();
    physical_data.set_window_size(win_width * 0.99, win_width * 0.99);
    physical_data.set_modes(2,-1);
    for (auto& floe : m_problem.get_floe_group().get_floes()){
        floe.static_floe().set_thickness(floe.static_floe().thickness() * max_size / 250);
    }
    real_type end_time = 500;
    bool init = true;
    do {
        m_problem.solve(end_time, 10, 10, init);
        m_problem.get_floe_group().stop_floes_in_window(win_width, win_width);
        std::cout << " Concentration : " << m_problem.floe_concentration() << std::endl;
        end_time += 500; init = false;
        if (*m_problem.QUIT) break; // exit normally after SIGINT
    } while (m_problem.get_floe_group().kinetic_energy() != 0 && end_time < 1e6);
    m_problem.get_floe_group().reset_impulses();
}

template<typename TProblem>
void
Generator<TProblem>::random_floe_group(std::size_t n, real_type max_size)
{
    auto& list_floes = m_problem.get_floe_group().get_floes();
    // auto sizes = random_size_repartition(n, max_size);
    auto sizes = exp_size_repartition(n, max_size);
    auto centers = spiral_distribution(sizes, max_size);
    auto generator = floe::random::get_uniquely_seeded_generator();
    std::uniform_int_distribution<int> distribution(0, m_biblio_size - 1);
    // unused pick to avoid first choice often being the same
    // (for quasi-simultaneous execs, even with different seeds)
    distribution(generator);
    list_floes.resize(n);

    for (std::size_t i = 0; i < n; i++)
    {
        int idx = distribution(generator);
        auto& base_shape = m_biblio_floe_h[idx];
        auto mesh = m_biblio_floe_h_meshes[idx];
        polygon_type shape;
        geometry::transform( base_shape, shape, scale_transformer<real_type>{ sizes[i] } );
        geometry::transform( mesh, mesh, scale_transformer<real_type>{ sizes[i] } );

        // Create Kinematic floe
        auto& floe = list_floes[i];
        // link static floe
        floe.attach_static_floe_ptr(std::unique_ptr<static_floe_type>(new static_floe_type()));
        
        auto& static_floe = floe.static_floe();

        // Import boundary
        std::unique_ptr<typename floe_type::geometry_type> geometry(new typename floe_type::geometry_type(shape));
        static_floe.attach_geometry_ptr(std::move(geometry));
        
        // Import mesh
        mesh_type& floe_mesh = floe.get_floe_h().m_static_mesh;
        floe_mesh = mesh;
        floe.static_floe().attach_mesh_ptr(&floe_mesh);
        // Done.

        m_problem.get_floe_group().get_floe_group_h().add_floe(floe.get_floe_h());

        // Set space-time state
        typename floe_type::state_type state;
        state.pos = centers[i];
        floe.set_state( state );
    }

}


template<typename TProblem>
std::vector<typename Generator<TProblem>::real_type>
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
std::vector<typename Generator<TProblem>::real_type>
Generator<TProblem>::exp_size_repartition(std::size_t n, real_type R_max)
{
    std::vector<real_type> v;
    real_type alpha = 1.5;
    real_type R_min = 0; // no min (resize if min floe too small)
    int nb_floes_per_size = 1;
    for (std::size_t i = 1; i <= n/nb_floes_per_size; i++)
    {
        // real_type R =  exp( (- 1 /  alpha) * log(i) + log(R_max));
        real_type R =  std::max(R_max * exp( (- 1 /  alpha) * log(i) ), R_min);
        for (int j=0; j<nb_floes_per_size; ++j){
            v.push_back(R);
        }
    }
    std::srand ( unsigned ( std::time(0) ) ); // seed for not having pseudorandom
    std::random_shuffle ( v.begin() + 1, v.end() ); // first floe(biggest) stays first (-> initially at center)
    return v;
}

template<typename TProblem>
std::vector<typename Generator<TProblem>::point_type>
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
            d_theta = r * 1.2 / R_max, d_R = r * 1.2;
        else
            d_theta = r * 1.2 / R, d_R = (Rmax / M_PI + 1) * (r / R);

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
