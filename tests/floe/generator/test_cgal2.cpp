#include "../tests/catch.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
// #include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <iostream>
#include <vector>
#include <cmath>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;


TEST_CASE( "Test Cgal constrained mesh generator", "[generator]" )
{
  CDT cdt;
  double R = 10;
  std::vector<Vertex_handle> v;
  cdt.insert( Point{0,0} );
  for (int i=0; i<10; ++i)
  {
      double angle = 2 * i *  M_PI / 20;
      v.push_back( cdt.insert(Point{R * cos(angle), R * sin(angle)}) );
  }
  for (auto it = v.begin(); it + 1!=v.end(); ++it)
  {
    cdt.insert_constraint(*it, *(it+1));
  }
  cdt.insert_constraint(*(v.end() - 1), *v.begin());

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
  std::cout << "Number of faces: " << cdt.number_of_faces() << std::endl;

  std::cout << "Meshing the triangulation..." << std::endl;
  const int nb_cells = 1000;
  const double cell_size{sqrt(2 * M_PI * R * R / nb_cells)}; // 
  CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, cell_size));
  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
  std::cout << "Number of faces: " << cdt.number_of_faces() << std::endl;

  /*
  CGAL::Unique_hash_map<Vertex_handle,int> V;
  int inum = 0;
  for( auto vit= cdt.vertices_begin(); vit != cdt.vertices_end() ; ++vit) {
      std::cout << vit->point().x() << ", " << vit->point().y() << std::endl;
      V[vit] = inum++;
  }
  // Faces to triangles
  for( auto ib = cdt.finite_faces_begin(); ib != cdt.finite_faces_end(); ++ib) {
      std::cout << V[ib->vertex(0)] << ", " << V[ib->vertex(1)] << ", " << V[ib->vertex(2)] << std::endl;
  }
  */
}
