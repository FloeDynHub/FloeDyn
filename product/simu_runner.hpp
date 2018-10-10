/*!
 * \file floe/product/simu_runner.hpp
 * \brief Simulation runner
 * \author Quentin Jouet
 */

#ifndef PRODUCT_SIMU_RUNNER_HPP
#define PRODUCT_SIMU_RUNNER_HPP
#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include "../product/config/interrupt.hpp"
#include "../product/config/config.hpp"
#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace types;
using namespace std;

namespace product {

/* Function used to check that of 'for_what' is specified, then
   'required_option' is specified too. */
void option_dependency(const po::variables_map& vm,
                        const char* for_what, const char* required_option)
{
    if (vm.count(for_what) && !vm[for_what].defaulted())
        if (vm.count(required_option) == 0 || vm[required_option].defaulted())
            throw logic_error(std::string("Option '") + for_what 
                              + "' requires option '" + required_option + "'.");
}

void handle_exception(std::exception& e){
    std::cerr << "Error: " << e.what() << "\n";
}

class SimuRunner{
public:
    SimuRunner( int argc, char* argv[] ) : desc("Allowed options"){
        this->init_program_options(argc, argv);
        this->init_interruption();
    }

    virtual int run(){
        if (vm.count("help")) {  
            cout << desc << "\n";
            return 0;
        }
        if (!this->check_options()){
            return 1;
        }
        #ifdef _OPENMP
        // omp_set_num_threads(1);
        Eigen::initParallel();
        #endif

        bool generate_floes = false;
        if (input_file_name == "generator") generate_floes = true; 

        problem_type P(epsilon, OBL_status);
        P.QUIT = &QUIT;
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
            if (!vm.count("nbfloes") || !vm.count("concentration"))
            {
                cerr << "Error : Generator requires --nbfloes and --concentration options.\n";
                return 1;
            }
            int nb_floes = vm["nbfloes"].as<int>();
            double concentration = vm["concentration"].as<value_type>();
            if (concentration <= 0 || concentration >= 1)
            {
                cerr << "Error : Generator constraint : 0 < concentration < 1.\n";
                return 1;
            }
            generator_type G( alpha, nbfpersize );
            G.set_exit_signal(&QUIT); // clean interrupt
            G.generate_floe_set(nb_floes, concentration, max_size, force_modes, force_speeds);
            P.set_floe_group(G.get_floe_group());
            P.get_floe_group().set_mu_static(0.7);
            #ifdef PBC
            auto win = P.get_floe_group().get_initial_window();
            P.set_topology(topology_type(win[0], win[1], win[2], win[3]));
            #endif
            P.make_input_file();
        }

        std::cout << "read TOPAZ" << std::endl;
        P.load_matlab_topaz_data(matlab_topaz_filename);
        P.get_dynamics_manager().set_rand_speed_add(rand_speed_add);
        P.get_dynamics_manager().set_norm_rand_speed(rand_norm);
        P.get_dynamics_manager().get_external_forces().get_physical_data().set_modes(force_modes[0],force_modes[1]);
        P.get_dynamics_manager().get_external_forces().get_physical_data().set_speeds(force_speeds[0],force_speeds[1]);
        // P.get_dynamics_manager().get_external_forces().get_physical_data().set_storm_mode(); // for simu: with storm
        // P.get_dynamics_manager().get_external_forces().get_physical_data().set_modes(2,0);   // for simu: ?
        // P.get_dynamics_manager().get_external_forces().get_physical_data().set_modes(-1,4);  // for simu: floes against obstacle

        #ifdef MULTIOUTPUT
            P.get_out_manager().set_size(nb_floe_select);
        #endif
        #ifdef LCPSTATS
            P.get_lcp_manager().get_solver().set_max_storage_sol(max_storage[0]);
            P.get_lcp_manager().get_solver().set_max_storage_unsol(max_storage[1]);
        #endif


        if (vm.count("rectime"))
        {
            P.recover_states_from_file(vm["recfile"].as<string>(), vm["rectime"].as<value_type>());
        }

        std::cout << "SOLVE..." << std::endl;
        // cout.precision(17);
        // std::cout << P.get_floe_group().total_area();
        P.get_floe_group().randomize_floes_thickness(random_thickness_coeff);
        P.get_floe_group().randomize_floes_oceanic_skin_drag(0.01);
        P.solve(endtime, default_time_step, out_time_step);
        return 0;
    }

protected:
    //!< Declare the supported options with BOOST
    po::options_description desc;
    po::variables_map vm;

    //!< supported options
    string input_file_name;
    #ifdef MULTIOUTPUT
        std::size_t nb_floe_select; //!< the size of the floe selection for the multiple output files (it is required!)
    #endif
    #ifdef LCPSTATS
        std::vector<int> max_storage; //!< solved and unsolved lcp max number 
    #endif
    value_type              endtime;
    value_type              default_time_step       = 10;
    value_type              out_time_step           = 60;
    std::vector<int>        force_modes             = std::vector<int>(2,0);
    std::vector<value_type> force_speeds            = std::vector<value_type>(2,0);
    int                     OBL_status              = 0;
    value_type              epsilon                 = 0.4;
    value_type              random_thickness_coeff  = 0.01;
    string                  matlab_topaz_filename   = "io/inputs/DataTopaz01.mat";
    value_type              max_size                = 250;
    bool                    rand_speed_add          = 0;
    value_type              rand_norm               = 1e-7;
    value_type              alpha                   = 1.5;
    int                     nbfpersize              = 1;

    void init_program_options( int argc, char* argv[] ){
        desc.add_options()
        // First parameter describes option name/short name
        // The second is parameter to option
        // The third is description
        ("help,h", "print usage message")
        ("input,i", po::value(&input_file_name)->required(), "input file path")
        ("fext, z", po::value(&matlab_topaz_filename)->default_value(matlab_topaz_filename), "external forces input file")
        #ifdef MULTIOUTPUT
            ("nbsefloes", po::value<std::size_t>(&nb_floe_select)->required(), "the size of the floe selection for the multiple output files")
        #endif
        #ifdef LCPSTATS
            ("lcpstats,l", po::value< std::vector<int> >(&max_storage)->required()->multitoken(), "both solved and unsolved lcp max number (2 entries)")
        #endif    

        ("fmodes", po::value< std::vector<int> >(&force_modes)->multitoken(), "forces modes [air, water].\n"
            "Possibilities: \n\n"

            "   (By default) NO atmospheric and water currents:\n"
            "       air mode: 0      water mode: 0\n\n"

            "   From a atmospheric and water fields: (by default)\n"
            "       air mode: 1      water mode: 1\n\n"

            "   For a storm (as a vortex): \n"
            "       air mode: 5      water mode: 0\n\n"

            "   For the initial floe pack generation: \n"
            "       air mode: 2      water mode: 0\n\n"
            "   or  air mode: 0     water mode: 2\n\n"

            "   For the simulation of percution against an obstacle: \n"
            "       air mode: 4      water mode: 0\n"
            "   or  air mode: 0     water mode: 4\n\n")

        ("fspeeds", po::value< std::vector<value_type> >(&force_speeds)->multitoken(), "forces speeds [air, water] (m/s).\n"
            "Possibilities: \n\n"

            "   From a atmospheric and water fields: (by default)\n"
            "       (no need for speeds)\n"
            "       air speed: 0      water speed: 0\n\n"

            "   For a storm (as a vortex): \n"
            "       (no need for speeds)\n"
            "       air speed: 0      water speed: 0\n\n"

            "   For the initial floe pack generation: \n"
            "       air speed >=10    water speed: 0\n\n"
            "   or  air speed: 0  1<= water speed <=4\n\n"

            "   For the simulation of percution against an obstacle: \n"
            "       air speed: 10     water speed: 0\n"
            "   or  air speed: 0      water speed: 1\n\n")

        ("bustle", po::value<bool>(&rand_speed_add), "1 to active the additional random floe velocities.")
        ("nbustle", po::value<value_type>(&rand_norm), "norm of the additional random floe velocities.")

        ("tend,t", po::value(&endtime)->required(), "simulation duration (seconds)")
        ("step,s", po::value(&default_time_step)->default_value(default_time_step), "default time step")
        ("outstep,o", po::value(&out_time_step)->default_value(out_time_step), "output time step")
        ("obl", po::value(&OBL_status)->default_value(OBL_status), "OBL status (0 or 1)")
        ("rectime,r", po::value<value_type>(), "time to recover states from")
        ("recfile,f", po::value<string>(), "file name to recover states from")

        ("nbfloes,n", po::value<int>(), "generator : how many floes ?")
        ("concentration,c", po::value<value_type>(), "generator : floes concentration (between 0 and 1)")
        ("maxsize,m", po::value(&max_size)->default_value(
            max_size, std::to_string(max_size)), "generator : floe max size (radius)")
        ("alpha,a", po::value<value_type>(&alpha), "fractal dimension for the distribution power law.")
        ("nbfpersize", po::value<int>(&nbfpersize), "number of floes per size for the distribution power law.")

        ("epsilon,e", po::value(&epsilon)->default_value(
            epsilon, std::to_string(epsilon)), "collision restitution coeff")
        ("sigma", po::value(&random_thickness_coeff)->default_value(
            random_thickness_coeff, std::to_string(random_thickness_coeff)),
            "Normal distribution coeff (sigma) for random ice thickness variation around 1m")
        ;
        try {
            po::store(po::parse_command_line(argc, argv, desc), this->vm);
            po::notify(this->vm);
        } catch ( const std::exception& e ) {
            std::cerr << e.what() << std::endl;
        }
    }

    void init_interruption(){ // WORKING ?
        // Handling interruption signal
        struct sigaction sa;
        memset( &sa, 0, sizeof(sa) );
        sa.sa_handler = interruption::got_signal;
        sigfillset(&sa.sa_mask);
        sigaction(SIGINT,&sa,NULL);
        sigaction(SIGSEGV,&sa,NULL); 
    }

    bool check_options(){
        try {
            po::notify(vm);
            option_dependency(vm, "rectime", "recfile");
            return true;
        } catch (std::exception& e) {
            handle_exception(e);
            return false;
        }
    }
};

} // namespace floe::product

#endif // PRODUCT_SIMU_RUNNER_HP