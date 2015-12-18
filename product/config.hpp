#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"
#include "floe/variable/floe_group.hpp"

#include "floe/problem/problem.hpp"

// #include "floe/ope/proximity_detector.hpp"
// #include "floe/collision/matlab/detector.h"

#include "floe/ope/LCP_manager.hpp"
#include "floe/ope/LCP_solver.h"
#include "floe/ope/generator_LCP_solver.hpp"
#include "floe/ope/collision_manager.hpp"

#include "floe/ope/proximity_detector.hpp"
#include "floe/collision/matlab/detector.h"

#include "floe/ope/physical_data.hpp"
#include "floe/ope/explicit_physical_data.hpp"
#include "floe/ope/external_forces.hpp"
#include "floe/ope/dynamics_manager.h"
#include "floe/domain/domain.hpp"

#include "floe/generator/generator.h"

#ifdef PBC
#include "floe/problem/periodic_problem.hpp"
#include "floe/topology/toric_topology.hpp"
#include "floe/collision/matlab/periodic_detector.h"
#include "floe/ope/periodic_dynamics_manager.h"
#endif



namespace ff = floe::floes;

using namespace floe::problem;

using real = double;
using floe_type = ff::KinematicFloe<ff::StaticFloe<real>>;
using floe_group_type = floe::variable::FloeGroup<floe_type>;

using proximity_detector_type = floe::ope::ProximityDetector<
    floe::collision::matlab::MatlabDetector<floe_type>
>;

using solver_type = floe::ope::LCPSolver<real>;
using manager_h_type = floe::ope::LCPManager<solver_type>;
using collision_manager_type = floe::ope::CollisionManager<manager_h_type>;
using generator_solver_type = floe::ope::GeneratorLCPSolver<real>;
using generator_manager_h_type = floe::ope::LCPManager<generator_solver_type>;
using generator_collision_manager_type = floe::ope::CollisionManager<generator_manager_h_type>;

using point_type = typename floe_type::point_type;
using physical_data_type = floe::ope::PhysicalData<point_type>;
using external_forces_type = floe::ope::ExternalForces<floe_type, physical_data_type>;
using generator_physical_data_type = floe::ope::ExplicitPhysicalData<point_type>;
using generator_external_forces_type = floe::ope::ExternalForces<floe_type, generator_physical_data_type>;
using generator_dynamics_manager_type = floe::ope::DynamicsManager<generator_external_forces_type, floe_group_type>;

using domain_type = floe::domain::Domain<real>;

#ifdef PBC // Periodic boundary conditions types
using topology_type = floe::topology::ToricTopology<typename floe_type::point_type>;
using periodic_proximity_detector_type = floe::ope::ProximityDetector<
    floe::collision::matlab::PeriodicMatlabDetector<floe_type, topology_type>
>;
using dynamics_manager_type = floe::ope::PeriodicDynamicsManager<external_forces_type, floe_group_type, topology_type>;
using problem_type = PeriodicProblem<
    floe_group_type,
    periodic_proximity_detector_type,
    collision_manager_type,
    dynamics_manager_type,
    domain_type,
    topology_type
>;
#else // Free boundary conditions types
using dynamics_manager_type = floe::ope::DynamicsManager<external_forces_type, floe_group_type>;
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
