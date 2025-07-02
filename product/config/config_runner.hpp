#ifndef PRODUCT_CONFIG_CONFIG_RUNNER_HPP
#define PRODUCT_CONFIG_CONFIG_RUNNER_HPP

#include "../product/simu_runner.hpp"

#ifdef MPIRUN
#include "../product/mpi_simu_runner.hpp"
#ifdef PBC
#include "../product/mpi_pbc_simu_runner.hpp"
#endif
#endif

namespace types {

#ifdef MPIRUN
    #ifdef PBC
        using simulation_runner_type = product::MPIPBCSimuRunner;
    #else
        using simulation_runner_type = product::MPISimuRunner;
    #endif
#else
    using simulation_runner_type = product::SimuRunner;
#endif

}

#endif // PRODUCT_CONFIG_CONFIG_RUNNER_HPP