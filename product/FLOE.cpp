#include "../product/simu_runner.hpp"
#include "../product/config/config_runner.hpp"
#include <iostream>
#include <chrono>

using namespace std::chrono;

int main( int argc, char* argv[] )
{
	auto t_start = high_resolution_clock::now();
    simulation_runner_type simu(argc, argv);
    simu.run();
	auto t_end = high_resolution_clock::now();

	auto simu_time = t_end-t_start;
	hours hh = duration_cast<hours> (simu_time);
	minutes mm = duration_cast<minutes> (simu_time % hours(1));
	seconds ss = duration_cast<seconds> (simu_time % minutes(1));
	std::cout 	<< "Chrono : Total simulation: " << hh.count() << "h " 
				<< mm.count() << "m " << ss.count() << "s \n";
    return 0;
}
