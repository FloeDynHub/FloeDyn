#include "../tests/catch.hpp"
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <iostream>


TEST_CASE( "Test blas", "[blas]" ) {

    using namespace std;
    boost::numeric::ublas::compressed_matrix<bool> adjacency(7, 9);
    adjacency(3,4)=1;
    for ( auto it1 = adjacency.begin1(); it1 != adjacency.end1(); ++it1 )
    {
        // cout << it1.index1();
        for ( auto it2 = it1.begin(); it2 != it1.end(); ++it2 )
            cout << it1.index1() << "." << it2.index2() << " ";
    }

}