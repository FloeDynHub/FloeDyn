#ifndef PRODUCT_CONFIG_CONFIG_RUNNER_HPP
#define PRODUCT_CONFIG_CONFIG_RUNNER_HPP

#include "../product/simu_runner.hpp"

#ifdef MPIRUN
#include "../product/mpi_simu_runner.hpp"
#endif

namespace types {

#ifdef MPIRUN
    using simulation_runner_type = product::MPISimuRunner;
#else
    using simulation_runner_type = product::SimuRunner;
#endif

}

#endif // PRODUCT_CONFIG_CONFIG_RUNNER_HPP