/*!
 * \file floe/product/mpi_simu_runner.hpp
 * \brief Simulation runner
 * \author Quentin Jouet
 */

#ifndef PRODUCT_MPI_SIMU_RUNNER_HPP
#define PRODUCT_MPI_SIMU_RUNNER_HPP
#include <mpi.h>
#include <iostream>
#include "../product/simu_runner.hpp"

namespace product {


class MPISimuRunner : public SimuRunner {
public:
    MPISimuRunner( int argc, char* argv[] ) : SimuRunner(argc, argv) {}

    virtual int run() override {
        if (this->vm.count("help")) {  
            cout << this->desc << "\n";
            return 0;
        }
        if (!this->check_options()){
            return 1;
        }
        #ifdef _OPENMP
        // omp_set_num_threads(1);
        Eigen::initParallel();
        #endif

        MPI_Init(nullptr, nullptr);
        int rank;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );

        int return_value;
        if (rank==0){
            // I'm the MASTER process
            std::cout << "MASTER OK" << std::endl;
            master_problem_type P(epsilon, OBL_status);
            return_value = this->run_problem(P);
        } else {
            // I'm a WORKER process
            std::cout << "WORKER #" << rank << " OK" << std::endl;
            worker_problem_type P(epsilon, OBL_status);
            return_value = this->run_problem(P);
        }
        MPI_Finalize();
        return return_value;
    }

private:

    template<typename TProblem>
    int run_problem(TProblem& P){
        P.QUIT = &QUIT;
        bool generate_floes = false;
        if (!generate_floes){
            try {
                P.load_config(this->vm["input"].as<string>());
            }
            catch(std::exception& e)
            {
                handle_exception(e);
                return 1;
            }
        }
        else {
            // todo generator
        }

        std::cout << "read TOPAZ" << std::endl;
        P.load_matlab_topaz_data(this->vm["fext"].as<string>());
        P.get_dynamics_manager().get_external_forces().get_physical_data().set_storm_mode();
        // To get same forcing as generator :
        // P.get_dynamics_manager().get_external_forces().get_physical_data().set_modes(2,0);
        // auto w = P.get_floe_group().get_initial_window();
        // P.get_dynamics_manager().get_external_forces().get_physical_data().set_window_size(w[1] - w[0], w[3] - w[2]);

        if (this->vm.count("rectime"))
        {
            P.recover_states_from_file(this->vm["recfile"].as<string>(), this->vm["rectime"].as<value_type>());
        }

        std::cout << "SOLVE..." << std::endl;
        P.get_floe_group().randomize_floes_thickness(this->vm["sigma"].as<value_type>());
        P.get_floe_group().randomize_floes_oceanic_skin_drag(0.01);
        P.solve(this->vm["tend"].as<value_type>(), this->vm["step"].as<value_type>(), this->vm["outstep"].as<value_type>());
        return 0;
    }
};

} // namespace floe::product

#endif // PRODUCT_SIMU_RUNNER_HP