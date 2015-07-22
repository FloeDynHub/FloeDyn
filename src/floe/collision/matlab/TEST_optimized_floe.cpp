/*!
 * \file floe/collision/matlab/TEST_optimized_floe.cpp
 * \brief Test OptimizedFloe class.
 * \author Roland Denis
 */

#include <iostream>
#include <cmath>
#include <cstddef>
#include <boost/math/constants/constants.hpp>

#include "floe/geometry/geometry.hpp"
#include "floe/floes/static_floe.hpp"
#include "floe/collision/matlab/optimized_floe.hpp"

int main()
{
    using namespace std;
    
    const size_t N = 10000;

    using real = double;
    using Floe = floe::floes::StaticFloe<real>;
    using geometry_type = typename Floe::geometry_type;

    const real pi = boost::math::constants::pi<real>();

    geometry_type* geo = new geometry_type{};
    auto& border = geo->outer();
    const real dtheta = -2*pi / N;
    for (size_t i = 0; i < N; ++i)
        border.push_back( { std::cos(i*dtheta), std::sin(i*dtheta) } );

    Floe floe;
    floe.attach_geometry_ptr(geo);
    
    cout << "Discrete circle area = " << floe.area() << endl;

    using OptimFloe = floe::collision::matlab::OptimizedFloe<Floe>;
    OptimFloe optim_floe(floe);

    cout << "cdist = " << optim_floe.cdist() << endl;
    cout << "tau = " << optim_floe.tau() << endl;
    cout << "global disk = " << floe::geometry::dsv( optim_floe.global_disk() ) << endl;
    //cout << "local disks = " << floe::geometry::dsv( optim_floe.local_disks() ) << endl;
    cout << "local disks count = " << optim_floe.local_disks().size() << endl;

    for (size_t i = 0; i <= 100*N; ++i)
    {
        floe.frame() = { {static_cast<real>(i), 0.}, i*pi/4 };
        optim_floe.update();
        //cout << "global disk = " << floe::geometry::dsv( optim_floe.global_disk() ) << endl;
        //cout << "local disks = " << floe::geometry::dsv( optim_floe.local_disks() ) << endl;
    }
    cout << "global disk = " << floe::geometry::dsv( optim_floe.global_disk() ) << endl;


    return 0;
}
