#include "../tests/catch.hpp"
#include <iostream>
#include "../product/config/config_dynamics.hpp"
#include <ctime>
#include <chrono>


namespace ff = floe::floes;

TEST_CASE( "Test Dynamics Manager", "[ope]" ) {

    using namespace std;
    value_type time = 0;
    dynamics_manager_type M(time);

    floe_group_type F;
    std::string mat_file_name = "tests/floe/io/matlab/1floe_rot1.mat";
    // std::string mat_file_name = "tests/floe/io/matlab/r1day_set_up_250sm_sz_60_list_so_350_str.mat";
    F.load_matlab_config(mat_file_name);

    cout << "EC begin : " << F.kinetic_energy() << endl;
    cout << "State begin : " << F.get_floes()[0].state() << endl;

    auto t_start = chrono::high_resolution_clock::now();
    std::clock_t start = std::clock();
    // for (int i = 0; i < 100; i++){
    //     M.move_floes(F, 0.5);
    // }
    auto t_end = chrono::high_resolution_clock::now();
    cout << "Chrono : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;

    cout << "EC end : " << F.kinetic_energy() << endl;
    cout << "State end : " << F.get_floes()[0].state() << endl;

}