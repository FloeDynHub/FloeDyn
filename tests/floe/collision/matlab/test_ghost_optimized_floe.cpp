#include "../tests/catch.hpp"

#include <iostream>
#include <cmath>
#include <cstddef>
#include <boost/math/constants/constants.hpp>

#include "floe/geometry/geometry.hpp"
#include "floe/floes/static_floe.hpp"
#include "floe/collision/matlab/optimized_floe.hpp"
#include "floe/collision/matlab/ghost_optimized_floe.hpp"

TEST_CASE( "Test ghost optimized floe", "[collision]" )
{
    using namespace std;
    
    const size_t N = 10000;

    using real = double;
    using Floe = floe::floes::StaticFloe<real>;
    using geometry_type = typename Floe::geometry_type;

    const real pi = boost::math::constants::pi<real>();

    auto geo = std::unique_ptr<geometry_type>(new geometry_type{});
    auto& border = geo->outer();
    const real dtheta = -2*pi / N;
    for (size_t i = 0; i < N; ++i)
        border.push_back( { std::cos(i*dtheta), std::sin(i*dtheta) } );

    Floe floe;
    floe.attach_geometry_ptr(std::move(geo));
    
    // cout << "Discrete circle area = " << floe.area() << endl;

    using OptimFloe = floe::collision::matlab::OptimizedFloe<Floe>;
    OptimFloe optim_floe(floe);

    // cout << "cdist = " << optim_floe.cdist << endl;
    // cout << "tau = " << optim_floe.tau << endl;
    // cout << "global disk = " << floe::geometry::dsv( optim_floe.global_disk() ) << endl;
    // cout << "local disks = " << floe::geometry::dsv( optim_floe.local_disks() ) << endl;
    // cout << "local disks count = " << optim_floe.local_disks().size() << endl;


    // cout << "global disk = " << floe::geometry::dsv( optim_floe.global_disk() ) << endl;


    using GhostOptim = floe::collision::matlab::GhostOptimizedFloe<OptimFloe>;
    using point_type = typename GhostOptim::point_type;
    point_type translation = point_type{3.5, 7};
    GhostOptim ghost( optim_floe, translation );

    REQUIRE( equal_points(ghost.global_disk().center, optim_floe.global_disk().center + translation) );
    REQUIRE(ghost.local_disks().size() == optim_floe.local_disks().size());

    vector<bool> eq(optim_floe.local_disks().size());
    for (std::size_t i = 0; i != eq.size(); ++i)
        eq[i] = equal_points(ghost.local_disks()[i].center, optim_floe.local_disks()[i].center + translation);

    REQUIRE( std::all_of(eq.begin(), eq.end(), [](bool b){ return b; }));

    REQUIRE( ghost.local_points() == optim_floe.local_points());


}
