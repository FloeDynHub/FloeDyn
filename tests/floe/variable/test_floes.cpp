#include "../tests/catch.hpp"
#include <iostream>

#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"
#include "floe/variable/floe_group.hpp"
#include "floe/geometry/arithmetic/arithmetic.hpp"
#include "floe/geometry/arithmetic/point_operators.hpp"

/* TEST Floe center of mass
#include "floe/integration/integrate.hpp"
#include "floe/integration/gauss_legendre.hpp"
*/

namespace ff = floe::floes;

TEST_CASE( "Test Floes", "[variable]" ) {

    using namespace floe::variable;
    using namespace std;

    using floe_type = ff::KinematicFloe<ff::StaticFloe<double>>;
    using point_type = typename floe_type::point_type;
    using mesh_type = floe_type::mesh_type;

    // using floe_alg_type = Floe_alg;
    // using floe_h_type = Floe_h<mesh_type, floe_alg_type>;
    // using floe_group_h_type = FloeGroup_h<floe_h_type>;

    // Create empty Floe Group
    // FloeGroup<floe_type, floe_group_h_type> F;
    FloeGroup<floe_type> F;
    REQUIRE( F.get_floes().size() == 0 );

    // Import floes from Matlab configuration
    std::string mat_file_name = "tests/floe/io/matlab/r1day_set_up_250sm_sz_60_list_so_350_str.mat";
    F.load_matlab_config(mat_file_name);
    REQUIRE(F.get_floes().size() == 350);
    // is discrete group correctly instanciated ?
    REQUIRE(F.get_floe_group_h().m_list_floe_h.size() == 350);


    /* TEST Floe center of mass
    using namespace floe::integration;
    using floe::integration::integrate;
    using real = double;
    const auto strategy = RefGaussLegendre<real,2,2>();

    for (std::size_t i=0; i< 20 ; i++)
    {
        floe_type& F0 = F.get_floes()[i];
        point_type V{0,0};
        for (auto& p : F0.mesh().points())
            V += p;

        V = (1. / F0.mesh().points().size()) * V;
        // cout << "Center of Mass " << V << endl;
        cout << "DIFF " << (V - F0.state().pos);
        auto result = integrate( [] (real x, real y) { return point_type{x,y}; }, F0.mesh(), strategy );
        cout << "DIFF2 " << (result / F0.area() - F0.state().pos) << endl;

    }
    */



}

// 19 192 466 242 004