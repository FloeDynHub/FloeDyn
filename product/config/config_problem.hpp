#ifndef PRODUCT_CONFIG_CONFIG_PROBLEM_HPP
#define PRODUCT_CONFIG_CONFIG_PROBLEM_HPP

#include "../product/config/config_floes.hpp"
#include "../product/config/config_detector.hpp"
#include "../product/config/config_collision.hpp"
#include "../product/config/config_dynamics.hpp"
#include "floe/problem/problem.hpp"

#ifdef PBC
#include "floe/problem/periodic_problem.hpp"
#endif

#ifdef MPIRUN
#include "floe/problem/mpi_master_problem.hpp"
#include "floe/problem/mpi_worker_problem.hpp"
#endif

namespace types {

using namespace floe::problem;

#ifdef PBC // Periodic boundary conditions types
using periodic_problem_type = PeriodicProblem<
    floe_group_type,
    periodic_proximity_detector_type,
    collision_manager_type,
    dynamics_manager_type,
    domain_type,
    topology_type
>;
#endif // Free boundary conditions types
using free_problem_type = Problem<
    floe_group_type,
    proximity_detector_type,
    collision_manager_type,
    dynamics_manager_type,
    domain_type
>;

#ifdef PBC
using problem_type = periodic_problem_type;
#else
using problem_type = free_problem_type;
#endif

#ifdef MPIRUN
#ifdef PBC
using master_problem_type = MPIMasterProblem<
    Problem<
        floe_group_type,
        master_proximity_detector_type,
        collision_manager_type,
        dynamics_manager_type,
        domain_type
    >
>;
using worker_grid_problem_type = MPIWorkerProblem<
    PeriodicProblem<
        floe_group_type,
        proximity_detector_type,
        collision_manager_type,
        dynamics_manager_type,
        domain_type,
        topology_type
    >
>;
using worker_in_border_problem_type = MPIWorkerProblem<free_problem_type>;
using worker_out_border_problem_type = MPIWorkerProblem<periodic_problem_type>;
#else
using master_problem_type = MPIMasterProblem<
    Problem<
        floe_group_type,
        master_proximity_detector_type,
        collision_manager_type,
        dynamics_manager_type,
        domain_type
    >
>;
#endif
using worker_problem_type = MPIWorkerProblem<free_problem_type>; // must be declared anyway for build
#else // NO MPI
#endif

} // namespace types

#endif // PRODUCT_CONFIG_CONFIG_PROBLEM_HPP