#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"
#include "floe/variable/floe_group.hpp"

#include "floe/problem/periodic_problem.hpp"

#include "floe/ope/proximity_detector.hpp"
#include "floe/topology/toric_topology.hpp"
#include "floe/collision/matlab/periodic_detector.hpp"

#include "floe/ope/collision_manager.hpp"
#include "floe/ope/periodic_dynamics_manager.hpp"
#include "floe/domain/domain.hpp"

#include "floe/integration/integrate.hpp"
#include "floe/integration/gauss_legendre.hpp"


namespace ff = floe::floes;

using namespace floe::problem;

using real = double;
using floe_type = ff::KinematicFloe<ff::StaticFloe<real>>;
using floe_group_type = floe::variable::FloeGroup<floe_type>;
using topology_type = floe::topology::ToricTopology<typename floe_type::point_type>;

using proximity_detector_type = floe::ope::ProximityDetector<
    floe::collision::matlab::PeriodicMatlabDetector<floe_type, topology_type>
>;

using collision_manager_type = floe::ope::CollisionManager<floe_type>;
using dynamics_manager_type = floe::ope::PeriodicDynamicsManager<floe_group_type, topology_type>;
using domain_type = floe::domain::Domain<real>;

using problem_type = PeriodicProblem<
    floe_group_type,
    proximity_detector_type,
    collision_manager_type,
    dynamics_manager_type,
    domain_type,
    topology_type
>;
