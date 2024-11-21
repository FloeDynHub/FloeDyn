/*!
 * \file floe/product/mpi_simu_runner.hpp
 * \brief Simulation runner
 * \author Quentin Jouet
 */

#ifndef PRODUCT_MPI_PBC_SIMU_RUNNER_HPP
#define PRODUCT_MPI_PBC_SIMU_RUNNER_HPP
#include <mpi.h>
#include <iostream>
#include "../product/mpi_simu_runner.hpp"

namespace product {


class MPIPBCSimuRunner : public MPISimuRunner {
public:
    MPIPBCSimuRunner( int argc, char* argv[] ) : MPISimuRunner(argc, argv) {}

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
        int total_processes;
        MPI_Comm_size(MPI_COMM_WORLD, &total_processes);
        int N = 1; // dim grid
        // iter from 0
        while (4 * N * N + 1 < total_processes) N++;
        if (rank==0){
            // I'm the MASTER process
            std::cout << "MASTER OK" << std::endl;
            master_problem_type P(epsilon, OBL_status);
            return run_problem(P);
        } else if (rank <= N * N){ 
            // I'm a WORKER grid process
            std::cout << "WORKER grid #" << rank << " OK" << std::endl;
            worker_grid_problem_type P(epsilon, OBL_status);
            return run_problem(P);
        } else if (rank <= 4 * N * (N-1) + 1){
            // I'm a WORKER in-border process (x-border or y-border or internal cross-border)
            std::cout << "WORKER in-border #" << rank << " OK" << std::endl;
            worker_in_border_problem_type P(epsilon, OBL_status);
            return run_problem(P);
        } else {
            // I'm a WORKER out-border process (periodic external border)
            std::cout << "WORKER out-border #" << rank << " OK" << std::endl;
            worker_out_border_problem_type P(epsilon, OBL_status);
            return run_problem(P);
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
        P.get_dynamics_manager().set_rand_speed_add(rand_speed_add);
        P.get_dynamics_manager().set_norm_rand_speed(rand_norm);
        P.get_dynamics_manager().get_external_forces().get_physical_data().set_modes(force_modes[0],force_modes[1]);
        P.get_dynamics_manager().get_external_forces().get_physical_data().set_speeds(force_speeds[0],force_speeds[1]);
        
        #ifdef MULTIOUTPUT
            P.get_out_manager().set_size(nb_floe_select);
        #endif
        #ifdef LCPSTATS
            P.get_lcp_manager().get_solver().set_max_storage_sol(max_storage[0]);
            P.get_lcp_manager().get_solver().set_max_storage_unsol(max_storage[1]);
        #endif

        if (this->vm.count("rectime"))
        {
            P.recover_states_from_file(this->vm["recfile"].as<string>(), this->vm["rectime"].as<value_type>());
        }

        //!< To get same forcing as generator (used only if force mode equal 2):
        if ((force_modes[0]==2 && force_modes[1]==0) || (force_modes[0]==0 && force_modes[1]==2)) {
            auto wi = P.get_floe_group().get_initial_window();
            std::cout << "initial windows is: " << wi[0] << ", " << wi[1] << ", " << wi[2] << ", " << wi[3] << "\n";
            auto wc = P.get_floe_group().bounding_window(0);
            std::cout << "current windows is: " << wc[0] << ", " << wc[1] << ", " << wc[2] << ", " << wc[3] << "\n";
            
            std::cout << "the initial concentration is: " << P.get_floe_group().initial_concentration() << "\n";
            std::cout << "the current concentration is: " << P.get_floe_group().floe_concentration() << "\n";

            value_type desired_conc;
            if ( vm.count("concentration") ) {
                desired_conc = vm["concentration"].as<value_type>();
            }
            else {desired_conc = 0.7;}
            std::cout << "The desired concentration is: " << desired_conc << "\n";

            value_type des_area = P.get_floe_group().total_area()/desired_conc;
            decltype(wc) wd;
            wd[0] = -std::sqrt(des_area)/2; wd[1] = -wd[0]; wd[2] = wd[0]; wd[3] = wd[1];
            std::cout << "desired windows is: " << wd[0] << ", " << wd[1] << ", " << wd[2] << ", " << wd[3] << "\n";

            P.get_dynamics_manager().get_external_forces().get_physical_data().set_window_size(wd[1] - wd[0], wd[3] - wd[2]);
            P.get_floe_group().stop_floes_in_window(wd[1] - wd[0], wd[3] - wd[2]);
        }

        std::cout << "SOLVE..." << std::endl;
        P.get_floe_group().set_mu_static(mu_static);
        if (mu_static!=0.7) {std::cout << "Warning: the ice/ice static friction coefficient is fixed to: " << mu_static << std::endl;}
        if (epsilon!=0.4) {std::cout << "Warning: the restitution coefficient is fixed to: " << epsilon << std::endl;}
        if (!P.get_floe_group().h5_contains_floes_characs(input_file_name)) {
            P.get_floe_group().randomize_floes_thickness(this->vm["sigma"].as<value_type>());
            P.get_floe_group().randomize_floes_oceanic_skin_drag(0.01);
        }
        P.solve(this->vm["tend"].as<value_type>(), this->vm["step"].as<value_type>(), this->vm["outstep"].as<value_type>());
        return 0;
    }
};

} // namespace floe::product

#endif // PRODUCT_MPI_PBC_SIMU_RUNNER_HPP