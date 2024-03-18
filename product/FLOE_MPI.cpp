#include <mpi.h>
#include <iostream>
#include <string>
#include <cassert>
#include "../product/config/interrupt.hpp"
#include "../product/config/config.hpp"
#include <boost/program_options.hpp>
namespace po = boost::program_options; 
using namespace std;
using namespace types;

/*
Main FloeDyn simulation with MPI parallelization (no PBC)
*/

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

template<typename TProblem>
int run_problem(TProblem& P, po::variables_map& vm){
    P.QUIT = &QUIT;
    bool generate_floes = false;
    if (!generate_floes){
        try {
            P.load_config(vm["input"].as<string>());
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
    P.load_matlab_topaz_data(vm["topaz"].as<string>());


    if (vm.count("rectime"))
    {
        P.recover_states_from_file("io/outputs/" + vm["recfile"].as<string>(), vm["rectime"].as<value_type>());
    }

    std::cout << "SOLVE..." << std::endl;
    // cout.precision(17);
    // std::cout << P.get_floe_group().total_area();
    P.get_floe_group().randomize_floes_thickness(vm["sigma"].as<value_type>());
    P.get_floe_group().randomize_floes_oceanic_skin_drag(0.01);
    P.solve(vm["tend"].as<value_type>(), vm["step"].as<value_type>(), vm["outstep"].as<value_type>());
    // MPI::Finalize();
    return 0;
}

int main( int argc, char* argv[] )
{
    string input_file_name;
    value_type endtime;
    value_type default_time_step = 10;
    value_type out_time_step = 60;
    int OBL_status = 0;
    value_type epsilon = 0.4;
    value_type random_thickness_coeff = 0.01;
    string matlab_topaz_filename = "io/library/DataTopaz01.mat";
    
    po::options_description desc("Allowed options");
    desc.add_options()
    // First parameter describes option name/short name
    // The second is parameter to option
    // The third is description
    ("help,h", "print usage message") 
    ("input,i", po::value(&input_file_name)->required(), "input file path")
    ("topaz", po::value(&matlab_topaz_filename)->default_value(matlab_topaz_filename), "wind/ocean data")
    ("tend,t", po::value(&endtime)->required(), "simulation duration (seconds)")
    ("step,s", po::value(&default_time_step)->default_value(default_time_step), "default time step")
    ("outstep,o", po::value(&out_time_step)->default_value(out_time_step), "output time step")
    ("obl", po::value(&OBL_status)->default_value(OBL_status), "OBL status (0 or 1)")
    ("rectime,r", po::value<value_type>(), "time to recover states from")
    ("recfile,f", po::value<string>(), "file name to recover states from")
    ("nbfloes,n", po::value<int>(), "generator : how many floes ?")
    ("concentration,c", po::value<value_type>(), "generator : floes concentration (between 0 and 1)")
    ("epsilon,e", po::value(&epsilon)->default_value(
        epsilon, std::to_string(epsilon)), "collision restitution coeff")
    ("sigma", po::value(&random_thickness_coeff)->default_value(
        random_thickness_coeff, std::to_string(random_thickness_coeff)),
        "Normal distribution coeff (sigma) for random ice thickness variation around 1m")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {  
        cout << desc << "\n";
        return 0;
    }

    try {
        po::notify(vm);
        option_dependency(vm, "rectime", "recfile");
    } catch (std::exception& e) {
        handle_exception(e);
        return 1;
    }

    // Handling interruption signal
    struct sigaction sa; 
    memset( &sa, 0, sizeof(sa) );
    sa.sa_handler = interruption::got_signal;
    sigfillset(&sa.sa_mask);
    sigaction(SIGINT,&sa,NULL);
    sigaction(SIGSEGV,&sa,NULL); 

    #ifdef _OPENMP
    // omp_set_num_threads(1);
    Eigen::initParallel();
    #endif

    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    if (rank==0){
        // I'm the MASTER process
        std::cout << "MASTER OK" << std::endl;
        master_problem_type P(epsilon, OBL_status);
        return run_problem(P, vm);
    } else {
        // I'm a WORKER process
        std::cout << "WORKER #" << rank << " OK" << std::endl;
        worker_problem_type P(epsilon, OBL_status);
        return run_problem(P, vm);
    }

}