#include "../tests/catch.hpp"
#include <iostream>

#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"
#include "floe/floes/floe_group.hpp"
#include "floe/collision/matlab/periodic_detector.hpp"
#include "floe/topology/toric_topology.hpp"
#include <cmath>
#include <chrono>

TEST_CASE( "Test Interpenetration", "[geometry]" ) {
    
    using namespace floe::floes;
    using namespace std;

    using value_type = double;
    using floe_type = floe::floes::KinematicFloe<floe::floes::StaticFloe<value_type>>;
    using point_type = typename floe_type::point_type;
    using TSpaceTopology = floe::topology::ToricTopology<point_type>;
    using TDetector = floe::collision::matlab::PeriodicMatlabDetector<floe_type, TSpaceTopology>;

    // Create floe group
    FloeGroup<floe_type> F;

    // Import floes from Matlab configuration
    std::string mat_file_name = "tests/floe/io/matlab/r1day_set_up_250sm_sz_60_list_so_350_str.mat";
    F.load_matlab_config(mat_file_name);

    // create topology (auto_topology() method in class PeriodicProblem)
    value_type min_x, min_y, max_x, max_y;
    min_x = min_y = std::numeric_limits<value_type>::max();
    max_x = max_y = - std::numeric_limits<value_type>::max();
    for (auto const& floe : F.get_floes())  
        for (auto const& pt : floe.geometry().outer())
        {
            min_x = std::min(min_x, pt.x);
            min_y = std::min(min_y, pt.y);
            max_x = std::max(max_x, pt.x);
            max_y = std::max(max_y, pt.y);
        }
    value_type margin = 1;
    auto space_topology = TSpaceTopology{min_x - margin, max_x + margin, min_y - margin, max_y + margin};


    // Create detector and link it to floes
    TDetector detector;
    detector.set_topology(space_topology);
    for (auto& floe : F.get_floes())
    {
        detector.push_back(&floe);
    }

    // update detector and do interpenetration test
    detector.update();
    bool I = detector.check_interpenetration();
    REQUIRE( I == true );

    // Recover states from an out file and do interpenetration test again
    /*
    F.recover_states_from_file("io/out15k.h5", 30293.7);
    detector.update();
    cout << "interpenetration test after recover : " << detector.check_interpenetration();
    */

}
