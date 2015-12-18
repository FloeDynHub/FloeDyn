#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"
#include "floe/floes/floe_group.hpp"

#include "floe/problem/problem.hpp"

// #include "floe/dynamics/proximity_detector.hpp"
// #include "floe/collision/matlab/detector.h"

#include "floe/lcp/LCP_manager.hpp"
#include "floe/lcp/solver/LCP_solver.h"
#include "floe/lcp/solver/generator_LCP_solver.hpp"
#include "floe/collision/collision_manager.hpp"

#include "floe/dynamics/proximity_detector.hpp"
#include "floe/collision/matlab/detector.h"

#include "floe/dynamics/physical_data.hpp"
#include "floe/dynamics/explicit_physical_data.hpp"
#include "floe/dynamics/external_forces.hpp"
#include "floe/dynamics/dynamics_manager.h"
#include "floe/domain/domain.hpp"

#include "floe/generator/generator.h"

#ifdef PBC
#include "floe/problem/periodic_problem.hpp"
#include "floe/topology/toric_topology.hpp"
#include "floe/collision/matlab/periodic_detector.h"
#include "floe/dynamics/periodic_dynamics_manager.h"
#endif



namespace ff = floe::floes;

using namespace floe::problem;

using real = double;
using floe_type = ff::KinematicFloe<ff::StaticFloe<real>>;
using floe_group_type = floe::floes::FloeGroup<floe_type>;

using proximity_detector_type = floe::dynamics::ProximityDetector<
    floe::collision::matlab::MatlabDetector<floe_type>
>;

using solver_type = floe::lcp::solver::LCPSolver<real>;
using manager_h_type = floe::lcp::LCPManager<solver_type>;
using collision_manager_type = floe::collision::CollisionManager<manager_h_type>;
using generator_solver_type = floe::lcp::solver::GeneratorLCPSolver<real>;
using generator_manager_h_type = floe::lcp::LCPManager<generator_solver_type>;
using generator_collision_manager_type = floe::collision::CollisionManager<generator_manager_h_type>;

using point_type = typename floe_type::point_type;
using physical_data_type = floe::dynamics::PhysicalData<point_type>;
using external_forces_type = floe::dynamics::ExternalForces<floe_type, physical_data_type>;
using generator_physical_data_type = floe::dynamics::ExplicitPhysicalData<point_type>;
using generator_external_forces_type = floe::dynamics::ExternalForces<floe_type, generator_physical_data_type>;
using generator_dynamics_manager_type = floe::dynamics::DynamicsManager<generator_external_forces_type, floe_group_type>;

using domain_type = floe::domain::Domain<real>;

#ifdef PBC // Periodic boundary conditions types
using topology_type = floe::topology::ToricTopology<typename floe_type::point_type>;
using periodic_proximity_detector_type = floe::dynamics::ProximityDetector<
    floe::collision::matlab::PeriodicMatlabDetector<floe_type, topology_type>
>;
using dynamics_manager_type = floe::dynamics::PeriodicDynamicsManager<external_forces_type, floe_group_type, topology_type>;
using problem_type = PeriodicProblem<
    floe_group_type,
    periodic_proximity_detector_type,
    collision_manager_type,
    dynamics_manager_type,
    domain_type,
    topology_type
>;
#else // Free boundary conditions types
using dynamics_manager_type = floe::dynamics::DynamicsManager<external_forces_type, floe_group_type>;
using problem_type = Problem<
    floe_group_type,
    proximity_detector_type,
    collision_manager_type,
    dynamics_manager_type,
    domain_type
>;
#endif

using generator_problem_type = Problem<
    floe_group_type,
    proximity_detector_type,
    generator_collision_manager_type,
    generator_dynamics_manager_type,
    domain_type
>;
using generator_type = floe::generator::Generator<generator_problem_type>;
