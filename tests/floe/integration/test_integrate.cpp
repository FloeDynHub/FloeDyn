/*!
 * \file floe/integration/TEST_integrate.cpp
 * \brief Test file for integrate overs a various set of geometrical objects.
 * \author Roland Denis
 */
#include "../tests/catch.hpp"
#include <iostream>

#include <cmath>

#include "floe/geometry/arithmetic/point_operators.hpp"

#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/adapted/c_array.hpp"
#include "floe/geometry/geometries/triangle.hpp"
#include "floe/geometry/geometries/point.hpp"
#include "floe/geometry/geometries/triangle_mesh.hpp"
#include "floe/integration/integrate.hpp"
#include "floe/integration/gauss_legendre.hpp"
#include "floe/geometry/core/cells.hpp"

BOOST_GEOMETRY_REGISTER_C_ARRAY_CS(cs::cartesian)

using namespace std;
using namespace floe;
using std::abs;

double cst_fun (double /* x */, double /* y */) {
    return 1.;
}

TEST_CASE( "Test Integrate", "[integration]" ) {
    using real = double;

    real eps = 1e-15;
    bool success;

    using namespace floe::integration;

    const auto strategy = RefGaussLegendre<real,2,2>();
    using floe::integration::integrate;

    {
    real result;
    const real triangle[3][2] = { {1.,0.}, {0.,1.}, {0.,0.} };

    //result = integrator.eval( [] (real x, real y) { return 1.; } );
    result = integrate( cst_fun, triangle, strategy );
    cout << "\\int_T 1 dx dy = " << result << " : " << ( ( success = abs(result - 0.5) <= eps ) ? "OK" : "KO" ) << endl;
    REQUIRE(success);
    
    result = integrate( [] (real x, real /* y */) { return x; }, triangle, strategy );
    cout << "\\int_T x dx dy = " << result << " : " << ( ( success = abs(result - 1./6.) <= eps ) ? "OK" : "KO" ) << endl;
    REQUIRE(success);

    result = integrate( [] (real /* x */ , real y) { return y; }, triangle , strategy);
    cout << "\\int_T y dx dy = " << result << " : " << ( ( success = abs(result - 1./6.) <= eps ) ? "OK" : "KO" ) << endl;
    REQUIRE(success);
    }

    {
    real result;
    const real triangle[3][2] = { {3.,0.}, {0.,1.}, {1.,0.} };

    result = integrate( cst_fun, triangle, strategy );
    cout << "\\int_T 1 dx dy = " << result << " : " << ( ( success = abs(result - 1.) <= eps ) ? "OK" : "KO" ) << endl;
    REQUIRE(success);

    result = integrate( [] (real x, real /* y */ ) { return x; }, triangle, strategy );
    cout << "\\int_T x dx dy = " << result << " : " << ( ( success = abs(result - 4./3.) <= eps ) ? "OK" : "KO" ) << endl;
    REQUIRE(success);

    result = integrate( [] (real /* x */, real y) { return y; }, triangle, strategy );
    cout << "\\int_T y dx dy = " << result << " : " << ( ( success = abs(result - 1./3.) <= eps ) ? "OK" : "KO" ) << endl;
    REQUIRE(success);

 
    }
    
    {
    real result;
    using Point = floe::geometry::Point<real>;
    using Triangle = floe::geometry::Triangle<Point>;
    const Triangle triangle = { {3.,0.}, {0.,1.}, {1.,0.} };

    result = integrate( cst_fun, triangle, strategy );
    cout << "\\int_T 1 dx dy = " << result << " : " << ( ( success = abs(result - 1.) <= eps ) ? "OK" : "KO" ) << endl;
    REQUIRE(success);

    result = integrate( [] (real x, real /* y */ ) { return x; }, triangle, strategy );
    cout << "\\int_T x dx dy = " << result << " : " << ( ( success = abs(result - 4./3.) <= eps ) ? "OK" : "KO" ) << endl;
    REQUIRE(success);

    result = integrate( [] (real /* x */, real y) { return y; }, triangle, strategy );
    cout << "\\int_T y dx dy = " << result << " : " << ( ( success = abs(result - 1./3.) <= eps ) ? "OK" : "KO" ) << endl;
    REQUIRE(success);

    Point pt_result = integrate( [] (real /* x */, real y) -> Point { return Point{y,2*y}; }, triangle, strategy );
    cout << "\\int_T (y,2y) dx dy = " << pt_result << " : " << ( ( success = norm2(pt_result - Point{1./3., 2./3.}) <= eps) ? "OK" : "KO" ) << endl;
    REQUIRE(success);
 
    }

    {
    real result;
    using Point = floe::geometry::Point<real>;
    using Mesh  = floe::geometry::TriangleMesh<Point>;

    Mesh mesh;

    mesh.points().push_back( { 0, 0 } );
    mesh.points().push_back( { 0, 1 } );
    mesh.points().push_back( { 1, 1 } );
    mesh.points().push_back( { 1, 0 } );
    
    mesh.add_triangle( 0, 2, 1 );
    mesh.add_triangle( 2, 0, 3 );

    result = integrate( [] (real /* x */, real /* y */) { return 1.; }, floe::geometry::cells(mesh)[0], strategy );
    cout << "\\int_M 1 dx dy = " << result << " : " << ( ( success = abs(result - 1./2.) <= eps ) ? "OK" : "KO" ) << endl;
    
    result = integrate( [] (real /* x */, real /* y */) { return 1.; }, mesh, strategy );
    cout << "\\int_M 1 dx dy = " << result << " : " << ( ( success = abs(result - 1.) <= eps ) ? "OK" : "KO" ) << endl;
    
    result = integrate( [] (real /* x */, real /* y */) { return 1.; }, mesh, strategy );
    cout << "\\int_M 1 dx dy = " << result << " : " << ( ( success = abs(result - 1.) <= eps ) ? "OK" : "KO" ) << endl;

    Point pt_result = integrate( [] (real /* x */, real y) { return Point{y, 2*y}; }, mesh, strategy );
    cout << "\\int_M (y,2y) dx dy = " << pt_result << " : " << ( ( success = norm2(pt_result - Point{0.5, 1.}) <= eps ) ? "OK" : "KO" ) << endl;

    pt_result = integrate( [] (real /* x */, real y) { return Point{0, 0}; }, mesh, strategy );
    cout << "\\int_M (0,0) dx dy = " << pt_result << " : " << ( ( success = norm2(pt_result - Point{0, 0}) <= eps ) ? "OK" : "KO" ) << endl;

    namespace trans = boost::geometry::strategy::transform;
    Mesh mesh2;
    trans::rotate_transformer<boost::geometry::degree, double, 2, 2> rotate(45.);
    floe::geometry::transform(mesh, mesh2, rotate);

    pt_result = integrate( [] (real /* x */, real y) { return Point{y, 2*y}; }, mesh2, strategy );
    cout << "\\int_M(45°) (y, 2y) dx dy = " << pt_result << " : " << ( ( success = norm2(pt_result - Point{0,0}) <= eps ) ? "OK" : "KO" ) << endl;
    

    }
}