/*!
 * \file floe/product/mpi_simu_runner.hpp
 * \brief Simulation runner
 * \author Quentin Jouet
 */

#ifndef PRODUCT_MPI_SIMU_RUNNER_HPP
#define PRODUCT_MPI_SIMU_RUNNER_HPP
#include <mpi.h>
#include <iostream>
#include <cstdlib> // setenv
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
            P.get_dynamics_manager().get_external_forces().set_O_latitude(O_latitude);
            return_value = this->run_problem(P);
        } else {
            // I'm a WORKER process
            std::cout << "WORKER #" << rank << " OK" << std::endl;
            worker_problem_type P(epsilon, OBL_status);
            P.get_dynamics_manager().get_external_forces().set_O_latitude(O_latitude);
            return_value = this->run_problem(P);
        }
        MPI_Finalize();
        return return_value;
    }

private:

    template<typename TProblem>
    int run_problem(TProblem& P){
        P.QUIT = &QUIT;
        bool generate_floes = (input_file_name == "generator");
        if (!generate_floes){
            try {
                P.load_config(input_file_name);
            }
            catch(std::exception& e)
            {
                handle_exception(e);
                return 1;
            }
        }
        else {
            // Only the MASTER generates the pack (sequentially, on a throwaway problem_type) and writes
            // it to an input .h5; the path is broadcast so EVERY process load_config()s the same file —
            // i.e. every rank holds the full floe geometries, exactly like a normal input run. (The
            // generation itself isn't distributed yet; that's the next brick via InterProcessMessage.)
            if (!vm.count("fmodes")) force_modes = {2, 0};
            int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            std::string fname;
            if (rank == 0) {
                problem_type genP(epsilon, OBL_status);
                genP.QUIT = &QUIT;
                genP.get_dynamics_manager().get_external_forces().set_O_latitude(O_latitude);
                fname = this->run_generator(genP);
            }
            int len = (int)fname.size();
            MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD); // workers block here until the master has written
            if (len == 0) return 1; // generation failed on the master -> everyone aborts
            fname.resize(len);
            MPI_Bcast(&fname[0], len, MPI_CHAR, 0, MPI_COMM_WORLD);
            // The file is written+closed by the master before the broadcasts above release the workers, so
            // all ranks now load it concurrently, exactly like a normal multi-process input run.
            try {
                P.load_config(fname);
            }
            catch(std::exception& e)
            {
                handle_exception(e);
                return 1;
            }
            input_file_name = fname; // from now on refer to the generated file, not the "generator" keyword
                                     // (e.g. the h5_contains_floes_characs check before solve opens it)
        }

        // Forcing: mirror the sequential runner — NetCDF data for mode 9/9, Matlab/TOPAZ for mode 1,
        // nothing otherwise (the old code loaded TOPAZ unconditionally from a non-existent "fext" option).
        if (force_modes[0] == 9 && force_modes[1] == 9) {
            P.get_dynamics_manager().get_external_forces().get_physical_data().load_nc_forcing_data(forcing_file_name);
        } else if (force_modes[0] == 1 || force_modes[1] == 1) {
            P.load_matlab_topaz_data(forcing_file_name);
        }
        P.get_dynamics_manager().set_rand_speed_add(rand_speed_add);
        P.get_dynamics_manager().set_norm_rand_speed(rand_norm);
        P.get_dynamics_manager().get_external_forces().get_physical_data().set_modes(force_modes[0],force_modes[1]);
        P.get_dynamics_manager().get_external_forces().get_physical_data().set_speeds(force_speeds[0],force_speeds[1]);
        // Honour --output in MPI (it was ignored, so the master kept its random default name). Only the
        // master writes output, so set the name on rank 0 only; workers never write and don't need it.
        {
            int out_rank; MPI_Comm_rank(MPI_COMM_WORLD, &out_rank);
            if (out_rank == 0 && !output_file_name.empty())
                P.get_out_manager().set_out_file_name(output_file_name);
        }

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
            P.get_floe_group().randomize_floes_oceanic_skin_drag(random_oceanic_skin_drag_coeff);
        }
        P.solve(this->vm["tend"].as<value_type>(), this->vm["step"].as<value_type>(), this->vm["outstep"].as<value_type>());
        return 0;
    }
};

} // namespace floe::product

#endif // PRODUCT_SIMU_RUNNER_HP