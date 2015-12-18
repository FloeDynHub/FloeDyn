#include "../tests/catch.hpp"

#include <iostream>

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
    read_list_so_from_file( mat_file_name, list_so);

    std::vector<TKinematicFloe> list_floes;
    list_so_to_floes( list_so, list_floes );

    REQUIRE(list_floes.size() == 350);

    for (auto& floe : list_floes)
        floe.update(); // TODO Avoid that ! (move doesn't move internal pointers)

    cout << list_floes[0].area() << endl;
    cout << fg::dsv(list_floes[0].geometry().outer()[0]) << endl;

}
