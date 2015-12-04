#include "../tests/catch.hpp"
#include <iostream>

// CGAL Mesh Generation
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

// CGAL Polygon approximation
#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>

// Boost geometry
#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/point.hpp"   
#include "floe/geometry/geometries/multi_point.hpp"

using point_type = floe::geometry::Point<double>;
using multi_point_type = floe::geometry::MultiPoint<point_type>;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_2<K>         Triangulation;
typedef Triangulation::Vertex_circulator Vertex_circulator;
typedef Triangulation::Point             Point;

TEST_CASE( "Test Cgal mesh generator", "[generator]" ) {

  Triangulation t;
  for (int i=0; i<8; ++i)
  {
      double angle = 2 * i *  M_PI / 8;
      t.insert(Point{cos(angle), sin(angle)});
  }

  // std::cout << t;
  std::cout << t.number_of_vertices () << std::endl;
  std::cout << t.number_of_faces () << std::endl;

  // for (auto i = t.finite_faces_begin (); i < t.finite_faces_end(); i++)
  //   std::cout << i->vertex(0)->point() << " " <<  i->vertex(1)->point() << " " << i->vertex(2)->point() << std::endl;

  // for (auto face : t.tds().faces())
  //   std::cout << face.vertex(0)->point() << " " <<  face.vertex(1)->point() << " " << face.vertex(2)->point() << std::endl;

  // for (auto i = t.finite_vertices_begin (); i < t.finite_vertices_end(); i++)
  //   std::cout << i->point() << std::endl;

  // for (auto i = t.all_edges_begin (); i < t.all_edges_end(); i++)
  //   std::cout << *i << "*" << std::endl;

  CGAL::Unique_hash_map<Triangulation::Vertex_handle,int> V;
  CGAL::Unique_hash_map<Triangulation::Face_handle,int> F;
 
  int inum = 0;
  
  // other vertices
  for( auto vit= t.vertices_begin(); vit != t.vertices_end() ; ++vit) {
    V[vit] = inum++;
    std::cout << vit->point().x() << std::endl;
  }

  // vertices of the faces
  inum = 0;
  int dim = (t.dimension() == -1 ? 1 :  t.dimension() + 1);
  for( auto ib = t.finite_faces_begin(); ib != t.finite_faces_end(); ++ib) {
    F[ib] = inum++;
    for(int j = 0; j < dim ; ++j) {
      std::cout << V[ib->vertex(j)] << " ";
    }
    std::cout << *ib ;
    std::cout << std::endl;
  }


  //// SIMPLIFY POLYGON :

  std::cout << "------------------" << std::endl;
  namespace PS = CGAL::Polyline_simplification_2;
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Polygon_2<K>                   Polygon_2;
  typedef PS::Stop_below_count_ratio_threshold Stop;
  typedef PS::Squared_distance_cost            Cost;
  Polygon_2 polygon;
  for (int i=0; i<8; ++i)
  {
    double angle = 2 * i *  M_PI / 8;
    polygon.push_back(Point{cos(angle), sin(angle)});
  }
  Cost cost;
  polygon = PS::simplify(polygon, cost, Stop(0.5));

  std::cout.precision(12);
  std::cout << polygon << std::endl;

}



