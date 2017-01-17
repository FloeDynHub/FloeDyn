/*!
 * \file generator/mesh_generator.hpp
 * \brief Mesh generator
 * \author Quentin Jouet
 */

#ifndef GENERATOR_MESH_GENERATOR_HPP
#define GENERATOR_MESH_GENERATOR_HPP

// CGAL Mesh Generation
#include <fstream>
#include <list>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>


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

template<typename TShape, typename TMesh>
TMesh generate_mesh_for_shape(TShape const& shape)
{
    using point_type = typename TShape::point_type;
    using value_type = typename point_type::value_type;
    CDT cdt;
    std::vector<Vertex_handle> v;
    for (point_type pt : shape.outer())
    {
      v.push_back( cdt.insert(Point{pt.x, pt.y}) );
    }
    for (auto it = v.begin(); it + 1!=v.end(); ++it)
    {
        cdt.insert_constraint(*it, *(it+1));
    }
    if (*(v.end() - 1) != *v.begin())
    {
        cdt.insert_constraint(*(v.end() - 1), *v.begin());
    }

    std::list<Point> list_of_seeds;
    list_of_seeds.push_back(Point(1000, 1000)); // /!\ domain must contain point (0,0)

    const int nb_cells = 20;
    const value_type cell_size{ sqrt(2 * floe::geometry::area(shape) / nb_cells) };
    // CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, cell_size));

    CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(), Criteria(0.125, cell_size));
    // std::cout << "Number of finite faces: " << cdt.number_of_faces() << std::endl;
    // int mesh_faces_counter = 0;
    // for(CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
    //   fit != cdt.finite_faces_end(); ++fit) 
    // {
    // if(fit->is_in_domain()) ++mesh_faces_counter;
    // }
    // std::cout << "Number of faces in the mesh domain: " << mesh_faces_counter << std::endl;


    // Convert to mesh_type
    TMesh mesh;
    CGAL::Unique_hash_map<Vertex_handle,int> V;

    int inum = 0;

    // Vertices to points
    auto& points = mesh.points();
    for( auto vit= cdt.vertices_begin(); vit != cdt.vertices_end() ; ++vit) {
        points.push_back(point_type{vit->point().x(), vit->point().y()});
        V[vit] = inum++;
    }
    // Faces to triangles
    for( auto ib = cdt.finite_faces_begin(); ib != cdt.finite_faces_end(); ++ib) {
        if(ib->is_in_domain())
            mesh.add_triangle(V[ib->vertex(0)], V[ib->vertex(1)], V[ib->vertex(2)]);
    }

    return mesh;
}


}} // namespace floe::generator


#endif // GENERATOR_MESH_GENERATOR_HPP
