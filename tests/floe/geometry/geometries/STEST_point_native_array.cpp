#include <iostream>

#include "geometry/primitive/point_native_array.hpp"

using namespace std;

int main() {
    using real = double;

    using geometry::primitive::Point;

    Point<real,3> pt;
    pt[0] = 1.;
    pt[1] = 2.;
    pt[2] = 3.;

    cout << pt[1] << endl;
    cout << get<1>(pt) << endl;

    return 0;
}
