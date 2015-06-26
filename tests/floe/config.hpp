
#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"
#include "floe/variable/floe_group.hpp"
#include "floe/ope/proximity_detector.hpp"
#include "floe/collision/matlab/detector.hpp"
#include "floe/problem/problem.hpp"
#include "floe/ope/collision_manager.hpp"
#include "floe/ope/dynamics_manager.hpp"
#include "floe/domain/domain.hpp"

#include "floe/integration/integrate.hpp"
#include "floe/integration/gauss_legendre.hpp"


namespace ff = floe::floes;

using namespace floe::problem;

using real = double;
using floe_type = ff::KinematicFloe<ff::StaticFloe<real>>;
using floe_group_type = floe::variable::FloeGroup<floe_type>;

using proximity_detector_type = floe::ope::ProximityDetector<
    floe::collision::matlab::MatlabDetector<floe_type>
>;

using collision_manager_type = floe::ope::CollisionManager;
using dynamics_manager_type = floe::ope::DynamicsManager<floe_group_type>;
using domain_type = floe::domain::Domain;

using problem_type = Problem<
    floe_group_type,
    proximity_detector_type,
    collision_manager_type,
    dynamics_manager_type,
    domain_type
>;
