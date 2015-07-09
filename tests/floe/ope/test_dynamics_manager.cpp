#include "../tests/catch.hpp"
#include <iostream>
#include "floe/ope/dynamics_manager.hpp"
#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"
#include "floe/variable/floe_group.hpp"
// #include <boost/timer/timer.hpp>
#include <ctime>
#include <chrono>


namespace ff = floe::floes;

TEST_CASE( "Test Dynamics Manager", "[ope]" ) {

    using namespace floe::ope;
    using namespace std;
    using real = double;
    using floe_type = ff::KinematicFloe<ff::StaticFloe<real>>;
    using floe_group_type = floe::variable::FloeGroup<floe_type>;
    DynamicsManager<floe_group_type> M;

    floe_group_type F;
    std::string mat_file_name = "tests/floe/io/matlab/1floe_rot1.mat";
    // std::string mat_file_name = "tests/floe/io/matlab/r1day_set_up_250sm_sz_60_list_so_350_str.mat";
    F.load_matlab_config(mat_file_name);

    cout << "EC begin : " << F.kinetic_energy() << endl;
    cout << "State begin : " << F.get_floes()[0].state() << endl;

    auto t_start = chrono::high_resolution_clock::now();
    std::clock_t start = std::clock();
    for (int i = 0; i < 10; i++){
        M.move_floes(F, 0.5);
    }
    // double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    auto t_end = chrono::high_resolution_clock::now();
    // cout << "CPU time : "<< 1000.0 * duration << " ms" << endl;
    cout << "Chrono : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;

    cout << "EC end : " << F.kinetic_energy() << endl;
    cout << "State end : " << F.get_floes()[0].state() << endl;


}