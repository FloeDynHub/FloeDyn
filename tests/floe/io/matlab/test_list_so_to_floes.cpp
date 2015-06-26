#include "../tests/catch.hpp"

#include <iostream>
#include <boost/timer/timer.hpp>

#include "floe/geometry/geometry.hpp"

#include "floe/io/matlab/list_so_to_floes.hpp"
#include "floe/io/matlab/list_so_import.hpp"
#include "floe/io/matlab/list_so.hpp"

#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"

TEST_CASE( "Test list_so_to_floes", "[io/matlab]" )
{

    using namespace std;
    using namespace floe::io::matlab;
    namespace fg = floe::geometry;
    namespace ff = floe::floes;
    
    using real = double;
    using TStaticFloe = ff::StaticFloe<real>;
    using TKinematicFloe = ff::KinematicFloe<TStaticFloe>;

    std::string mat_file_name = "tests/floe/io/matlab/r1day_set_up_250sm_sz_60_list_so_350_str.mat";

    MatlabListSolid<double> list_so;
    cout << "Reading \"" << mat_file_name << "\" ... " << endl;
    boost::timer::cpu_timer timer;
    read_list_so_from_file( mat_file_name, list_so);
    timer.stop();
    cout << "\t" << timer.format() << endl;

    cout << "Importing floes ... " << endl;
    timer.start();
    auto floes = list_so_to_floes<TKinematicFloe>( list_so );
    timer.stop();
    cout << "\t" << timer.format() << endl;

    REQUIRE(floes.size() == 350);

    for (auto& floe : floes)
        floe.update(); // TODO Avoid that ! (move doesn't move internal pointers)

    cout << floes[0].area() << endl;
    cout << fg::dsv(floes[0].geometry().outer()[0]) << endl;

}
