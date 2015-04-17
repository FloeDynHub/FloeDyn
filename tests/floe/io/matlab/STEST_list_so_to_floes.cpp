/*!
 * \file floe/io/matlab/TEST_list_so_to_floes.cpp
 * \brief TEST file for list_so import from Matlab file.
 * \author Roland Denis
 */

#include <iostream>

#include <boost/timer/timer.hpp>

#include "floe/geometry/geometry.hpp"

#include "floe/io/matlab/list_so_to_floes.hpp"
#include "floe/io/matlab/list_so_import.hpp"
#include "floe/io/matlab/list_so.hpp"

#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"

int main( int argc, char* argv[] )
{
    using namespace std;
    using namespace floe::io::matlab;
    namespace fg = floe::geometry;
    namespace ff = floe::floes;
    
    using real = double;
    using TStaticFloe = ff::StaticFloe<real>;
    using TKinematicFloe = ff::KinematicFloe<TStaticFloe>;

    if ( argc < 2 )
    {
        cout << "Usage: " << argv[0] << " <matlab_file_name>" << endl;
        return 1;
    }

    MatlabListSolid<double> list_so;
    cout << "Reading \"" << argv[1] << "\" ... " << endl;
    boost::timer::cpu_timer timer;
    read_list_so_from_file( argv[1], list_so);
    timer.stop();
    cout << "\t" << timer.format() << endl;

    cout << "Importing floes ... " << endl;
    timer.start();
    auto floes = list_so_to_floes<TKinematicFloe>( list_so );
    timer.stop();
    cout << "\t" << timer.format() << endl;

    cout << floes.size() << endl;

    cout << floes[0]->area() << endl;
    cout << fg::dsv(floes[0]->geometry().outer()[0]) << endl;

    return 0;
}

