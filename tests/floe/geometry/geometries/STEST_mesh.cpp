#include <iostream>
#include <typeinfo>

#include <boost/concept_check.hpp>

#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/point.hpp"
#include "floe/geometry/geometries/triangle_mesh.hpp"

#include "floe/geometry/geometries/concepts/mesh_concept.hpp"

#include "floe/geometry/core/point_type.hpp"
#include "floe/geometry/core/coordinate_type.hpp"

int main() {
    using namespace std;
    namespace bg = boost::geometry;
    namespace fg = floe::geometry;

    using real  = double;
    using Point = fg::Point<real>;
    using Mesh  = fg::TriangleMesh<Point>;

    Mesh mesh;
    
    BOOST_CONCEPT_ASSERT( (fg::concepts::Mesh<Mesh>) );
    BOOST_CONCEPT_ASSERT( (fg::concepts::ConstMesh<Mesh>) );

    mesh.points().push_back( { 0, 0 } );
    mesh.points().push_back( { 0, 1 } );
    mesh.points().push_back( { 1, 1 } );
    mesh.points().push_back( { 1, 0 } );
    
    mesh.add_triangle( 0, 1, 2 );
    mesh.add_triangle( 2, 3, 0 );

    cout << fg::dsv( mesh.cells()[0] ) << endl;

    cout << fg::dsv( mesh.cells() ) << endl;

    typedef typename fg::point_type<Mesh>::type pt_type;
    cout << typeid(pt_type).name() << endl;

    typedef typename fg::coordinate_type<Mesh>::type c_type;
    cout << typeid(c_type).name() << endl;

    // Test transformations
    namespace trans = boost::geometry::strategy::transform;

    Mesh mesh2;
    trans::rotate_transformer<boost::geometry::degree, double, 2, 2> rotate(45.);
    fg::transform(mesh, mesh2, rotate);
    cout << fg::dsv( mesh2.cells() ) << endl;

    return 0;
}
