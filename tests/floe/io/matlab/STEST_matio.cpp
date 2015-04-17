#include <iostream>

#include "floe/io/matlab/list_so_import.hpp"

int main( int argc, char* argv[] )
{
    using namespace floe::io::matlab;
    MatlabListSolid<double> list_so;
    read_list_so_from_file(argv[1], list_so);

    return 0;
}
