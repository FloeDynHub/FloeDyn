/*!
 * \file floe/collision/matlab/TEST_detector.cpp
 * \brief Test file for the collision detector with random circular floes.
 * \author Roland Denis
 */

#include <iostream>

#include <boost/timer/timer.hpp>

#include "floe/geometry/geometry.hpp"
#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"
#include "floe/collision/matlab/detector.h"

int main(int argc, char* argv[])
{
    using namespace std;

    using real = double;
    using Floe = floe::floes::StaticFloe<real>;
    using MovingFloe = floe::floes::KinematicFloe<Floe>;
    using geometry_type = typename Floe::geometry_type;
    using Detector = floe::collision::matlab::MatlabDetector<MovingFloe>;
    
    const real pi = boost::math::constants::pi<real>();

    if (argc < 4)
    {
        cout << "Usage: " << argv[0] << " <N_floes> <N_point_per_floes> <domain_size>" << endl;
        return 1;
    }
    const size_t N_floes = atoi(argv[1]); // 2000;
    const size_t N_points = atoi(argv[2]); // 100
    const double L = atof(argv[3]); // 200;
    
    cout << "Parameters:" << endl;
    cout << "\t" << N_floes << " floes." << endl;
    cout << "\t" << N_points << " points per floe." << endl;
    cout << "\tDomain of size " << L << "x" << L << endl;
    cout << endl;

    // Creating floes
    std::vector<MovingFloe*> floe_list;
    cout << "Creating floe ... " << endl;
    boost::timer::cpu_timer timer;
    for (size_t i = 0; i < N_floes; ++i)
    {
        geometry_type* geo = new geometry_type{};
        auto& border = geo->outer();
        const real dtheta = -2*pi / N_points;
        for (size_t i = 0; i < N_points; ++i)
            border.push_back( { std::cos(i*dtheta), std::sin(i*dtheta) } );

        Floe* floe = new Floe();
        floe->attach_geometry_ptr(geo);

        MovingFloe* mfloe = new MovingFloe();
        mfloe->attach_static_floe_ptr(floe);

        floe_list.push_back(mfloe);

        //cout << mfloe->area() << endl;
    }
    timer.stop();
    cout << "\t" << timer.format() << endl;

    // Initializing detector
    cout << "Linking floes to detector ..." << endl;
    timer.start();
    Detector detector;
    for (auto& floe_ptr : floe_list)
        detector.push_back(floe_ptr);
    timer.stop();
    cout << "\t" << timer.format() << endl;


    // Moving floes
    cout << "Moving floes ..." << endl;
    timer.start();
    for ( std::size_t i = 0; i < floe_list.size(); ++i )
    {
        auto floe_ptr = floe_list[i];
        // floe_ptr->state().pos += {double(std::rand())*L/RAND_MAX,double(std::rand())*L/RAND_MAX};
        floe_ptr->state().pos = {2.*(i%80), 2.*(i/80)};
        floe_ptr->update();
    }
    timer.stop();
    cout << "\t" << timer.format() << endl;


    // Update collisions
    cout << "Updating optimizers and detecting contacts ..." << endl;
    timer.start();
    for (std::size_t i = 0; i < 1; ++i)
        detector.update();
    timer.stop();
    cout << "\t" << timer.format() << endl;

    cout << "Informations:" << endl;
    cout << "\t" << detector.num_points() << " points." << endl;
    cout << "\t" << detector.num_local_disks() << " local disks." << endl;
    // cout << "\t" << detector.num_contacts() << " contacts (" << (detector.num_contacts()*100./detector.num_points()) << "%)" << endl;
    cout << endl;

    // Freeing memory
    for ( auto& floe_ptr : floe_list )
        delete floe_ptr;
    

    return 0;

}
