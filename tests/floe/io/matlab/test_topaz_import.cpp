#include "../tests/catch.hpp"
#include <iostream>
#include "floe/io/matlab/topaz_import.hpp"
#include "floe/geometry/geometries/point.hpp" 


TEST_CASE( "Test Collision Manager", "[ope]" ) {

	using namespace floe::io::matlab;
	using point_vector_type = std::vector<floe::geometry::Point<double>>;

    std::string file_name = "tests/floe/io/matlab/DataTopaz01.mat";
    point_vector_type ocean_data;
    point_vector_type air_data;
    read_topaz_from_file(file_name, ocean_data, air_data);

    REQUIRE( ocean_data.size() == 696 );
    REQUIRE( air_data.size() == 696 );

}