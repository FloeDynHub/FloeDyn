/*!
 * \file floe/integration/TEST_gauss_legendre.cpp
 * \brief Test file for Gauss-Legendre quadrature method.
 * \author Roland Denis
 */

#include <iostream>
#include <cmath>

#include "gauss_legendre.hpp"

using namespace std;
using std::abs;

double cst_fun (double x, double y) {
    return 1.;
}

int main() {
    using real = double;

    real eps = 1e-15;
    real result;
    bool success;


    const auto integrator = integration::RefGaussLegendre<real,2,2>();

    //result = integrator.eval( [] (real x, real y) { return 1.; } );
    result = integrator.eval( cst_fun );
    cout << "\\int_T 1 dx dy = " << result << " : " << ( ( success = abs(result - 0.5) <= eps ) ? "OK" : "KO" ) << endl;
    if ( !success ) return 1;

    result = integrator.eval( [] (real x, real y) { return x; } );
    cout << "\\int_T x dx dy = " << result << " : " << ( ( success = abs(result - 1./6.) <= eps ) ? "OK" : "KO" ) << endl;
    if ( !success ) return 1;

    result = integrator.eval( [] (real x, real y) { return y; } );
    cout << "\\int_T y dx dy = " << result << " : " << ( ( success = abs(result - 1./6.) <= eps ) ? "OK" : "KO" ) << endl;
    if ( !success ) return 1;

    return 0;
}
