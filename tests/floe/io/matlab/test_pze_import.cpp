#include "../tests/catch.hpp"
#include <iostream>
#include "floe/geometry/geometries/point.hpp" 
#include "floe/io/matlab/pze_import.hpp"
#include "floe/topology/toric_topology.hpp"


TEST_CASE( "Test Collision Manager", "[ope]" ) {

	using namespace floe::io::matlab;
	using point_type = floe::geometry::Point<double>;
	using topology_type = floe::topology::ToricTopology<point_type>;

    std::string file_name = "io/Configs_09_2015/set_up_250sm_sz_30_list_so_350_disk_str.mat";

    topology_type T{read_pze_from_file<topology_type>(file_name)};
    
    std::cout << T.center();

}