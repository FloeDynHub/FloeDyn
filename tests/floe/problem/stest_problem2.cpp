#include <iostream>
#include "../tests/floe/config.hpp"


namespace floe { namespace problem{

template <
    typename TFloe,
    typename TProxymityDetector,
    typename TCollisionManager,
    typename TDynamicsManager,
    typename TDomain
>
void Problem<TFloe,TProxymityDetector,TCollisionManager,
             TDynamicsManager, TDomain>::test(){
}

}} // namespace floe::problem

int main( int argc, char* argv[] )
{
//    std::string mat_file_name = "/Users/Serge/Projects/FloeCpp/tests/floe/io/matlab/r1day_set_up_250sm_sz_60_list_so_350_str.mat";
    // std::string mat_file_name = "tests/floe/io/matlab/matlabv6.mat";
     std::string mat_file_name = "/Users/Serge/Projects/FloeCpp/tests/floe/io/matlab/config_q2_str.mat";

    problem_type P;
    P.load_matlab_config(mat_file_name);

    P.test();
    P.solve(10, 1);

    return 0;
}