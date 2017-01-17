#include "../tests/catch.hpp"
#include <iostream>
#include "floe/ope/time_scale_manager.hpp"
#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"
#include "floe/variable/floe_group.hpp"
#include "floe/collision/matlab/detector.hpp"
#include "floe/domain/domain.hpp"
#include <chrono>


namespace ff = floe::floes;

TEST_CASE( "Test Dynamics Manager", "[ope]" ) {

    using namespace std;
    using namespace floe::ope;
    using namespace std;
    using real = double;
    using floe_type = ff::KinematicFloe<ff::StaticFloe<real>>;
    using floe_group_type = floe::variable::FloeGroup<floe_type>;
    using TDetector = floe::collision::matlab::MatlabDetector<floe_type>;
    using TDomain = floe::domain::Domain;
    using time_scale_manager_type = TimeScaleManager<TDomain, TDetector>;

    TDomain domain;

    floe_group_type F;
    // std::string mat_file_name = "tests/floe/io/matlab/matlabv6.mat";
    std::string mat_file_name = "tests/floe/io/matlab/r1day_set_up_250sm_sz_60_list_so_350_str.mat";
    F.load_matlab_config(mat_file_name);

    TDetector detector;
    for (auto& floe : F.get_floes())
    {
        floe.update(); // TODO Avoid that ! (move doesn't move internal pointers)
        detector.push_back(&floe);
    }
    detector.update();

    // time_scale_manager_type M{domain, detector};
    time_scale_manager_type M;

    auto t_start = chrono::high_resolution_clock::now();
    cout << M.delta_t_secu(&domain, &detector) << endl;
    auto t_end = chrono::high_resolution_clock::now();
    cout << "Chrono : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;


}