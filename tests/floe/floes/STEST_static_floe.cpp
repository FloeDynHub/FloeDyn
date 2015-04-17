/*!
 * \file floe/floes/TEST_static_floe.cpp
 * \brief TEST file for StaticFloe and KinematicFloe classes.
 * \author Roland Denis
 */

#include <iostream>
#include <boost/geometry/geometries/polygon.hpp>
#include "floe/geometry/geometry.hpp"
#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"
#include "floe/floes/identifiable_mixin.hpp"
#include "floe/geometry/arithmetic/point_operators.hpp"

int main()
{
    namespace fg = floe::geometry;
    namespace ff = floe::floes;
    using real = double;
    
    using StaticFloe = ff::StaticFloe<real>;
    //using StaticFloe = ff::Identifiable<std::size_t, ff::StaticFloe<real,fg::Point<real>,boost::geometry::model::polygon<fg::Point<real> > > >;
    using Point = typename StaticFloe::point_type;

    StaticFloe* floe = new StaticFloe();

    typedef typename StaticFloe::geometry_type geometry_type;
    geometry_type* geo = new geometry_type();
    
    //auto& border = geo->boundary();
    auto& border = geo->outer();
    border.push_back({0.,0.});
    border.push_back({0.,1.});
    border.push_back({2.,2.});
    border.push_back({1.,0.});

    floe->attach_geometry_ptr(geo);

    using namespace std;
    cout << floe->area() << endl;

    using KinematicFloe = ff::KinematicFloe<StaticFloe>;
    KinematicFloe moving_floe;
    cout << moving_floe.area() << endl;

    moving_floe.attach_static_floe_ptr(floe);
    moving_floe.update();

    
    cout << fg::dsv(moving_floe.geometry()) << endl;
    moving_floe.state().pos = {1., 1.};
    moving_floe.update();
    cout << fg::dsv(moving_floe.geometry()) << endl;

    moving_floe.state().pos += Point{1.,0.};
    moving_floe.update();
    cout << fg::dsv(moving_floe.geometry()) << endl;
    return 0;
}
