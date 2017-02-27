#include "../product/simu_runner.hpp"
#include "../product/config/config_runner.hpp"

int main( int argc, char* argv[] )
{
    simulation_runner_type simu(argc, argv);
    return simu.run();
}