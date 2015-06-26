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
    using namespace std;

    if ( argc < 4 )
    {
        cout << "Usage: " << argv[0] << " <matlab_file_name> <end_time> <out_step_nb>" << endl;
        return 1;
    }

    std::string mat_file_name = argv[1];

    problem_type P;
    P.load_matlab_config(mat_file_name);

    P.solve(atoi(argv[2]), atof(argv[3]));

    return 0;
}