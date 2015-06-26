#include "../tests/catch.hpp"
#include <iostream>
#include "floe/problem/problem.hpp"
#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"
#include "floe/ope/proximity_detector.hpp"
#include "floe/collision/matlab/detector.hpp"

namespace ff = floe::floes;

TEST_CASE( "Test Problem", "[problem]" ) {
    using namespace floe::problem;

    using floe_type = ff::KinematicFloe<ff::StaticFloe<double>>;
    using proximity_detector_type = floe::ope::ProximityDetector<floe::collision::matlab::MatlabDetector<floe_type>>;
    using problem_type = Problem<floe_type, proximity_detector_type>;
    std::string mat_file_name = "tests/floe/io/matlab/r1day_set_up_250sm_sz_60_list_so_350_str.mat";

    problem_type P;
    P.load_matlab_config(mat_file_name);
    P.solve();
}