#include "../tests/catch.hpp"
#include <iostream>
#include "floe/geometry/geometries/point.hpp" 
#include "floe/io/matlab/pze_import.hpp"


TEST_CASE( "Test matlab pze import (ocean window)", "[io]" ) {

	using namespace floe::io::matlab;
	using value_type = double;
	using point_type = floe::geometry::Point<value_type>;

    std::string file_name = "io/inputs/Configs_09_2015/set_up_250sm_sz_30_list_so_350_disk_str.mat";

    auto T = read_pze_from_file(file_name);
    
    std::cout << T[0] << ", " << T[1] << ", " << T[2] << ", " << T[3];

}