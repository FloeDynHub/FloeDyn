#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"
#include "floe/ope/proximity_detector.hpp"
#include "floe/collision/matlab/detector.hpp"
#include "floe/problem/problem.hpp"
#include "floe/ope/collision_manager.hpp"


namespace ff = floe::floes;

using namespace floe::problem;

using real = double;
using floe_type = ff::KinematicFloe<ff::StaticFloe<double>>;

using proximity_detector_type = floe::ope::ProximityDetector<
    floe::collision::matlab::MatlabDetector<floe_type>
>;

using collision_manager_type = floe::ope::CollisionManager;

using problem_type = Problem<
    floe_type,
    proximity_detector_type,
    collision_manager_type
>;