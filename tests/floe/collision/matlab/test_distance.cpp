#include "../tests/catch.hpp"
#include <iostream>
#include <fstream>

#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"
#include "floe/variable/floe_group.hpp"
#include "floe/collision/matlab/periodic_detector.hpp"
#include "floe/topology/toric_topology.hpp"

// For debug, shouldn't be here
using VALUE_TYPE = double;
VALUE_TYPE DT_DEFAULT{100};
#include "floe/ope/time_scale_manager.hpp"
#include "floe/domain/domain.hpp"
// For debug, shouldn't be here

#include <cmath>
#include "matio.h"

using namespace std;


TEST_CASE( "Test distance between floes", "[collision]" ) {
    
    using namespace floe::variable;
    using namespace std;

    using value_type = double;
    using floe_type = floe::floes::KinematicFloe<floe::floes::StaticFloe<value_type>>;
    using point_type = typename floe_type::point_type;
    using TSpaceTopology = floe::topology::ToricTopology<point_type>;
    using TDetector = floe::collision::matlab::PeriodicMatlabDetector<floe_type, TSpaceTopology>;

    // Create floe group
    FloeGroup<floe_type> F;

    // Import floes from Matlab configuration
    // std::string mat_file_name = "io/set_up_250sm_sz_60_list_so_350_str.mat";
    std::string mat_file_name = "tests/floe/io/matlab/r1day_set_up_250sm_sz_60_list_so_350_str.mat";
    F.load_matlab_config(mat_file_name);

    // create topology 
    auto space_topology = TSpaceTopology{-600, 600, -600, 600}; // Pze data (in matlab config)

    // Create detector and link it to floes
    TDetector detector;
    detector.set_topology(space_topology);
    for (auto& floe : F.get_floes())
    {
        detector.push_back(&floe);
    }

    // update detector
    detector.update();

    
    auto const& d_secu = detector.get_dist_secu();
    auto const& d_opt = detector.get_dist_opt();
    /*
    // disp min d_secu
    value_type min_d = std::numeric_limits<value_type>::max();
    for (std::size_t n1 = 0; n1 < d_secu.size1(); ++n1)
        for (std::size_t n2 = n1+1; n2 < d_secu.size2(); ++n2)
            min_d = std::min(min_d, d_secu(n1, n2));
    cout << "MIN D_SECU = " << min_d << endl;
    
    // write d_secu in csv format
    ofstream myfile;
    myfile.open ("io/d_secu.csv");
    for (std::size_t n1 = 0; n1 < d_secu.size1(); ++n1)
    {
        for (std::size_t n2 = 0; n2 < d_secu.size2(); ++n2)
        {
            myfile << d_secu(n1, n2) << ",";
            min_d = std::min(min_d, d_secu(n1, n2));
        }
        myfile << "\n";
    }
    myfile.close();


    // write d_opt in csv format and count non-zero values
    myfile.open ("io/d_opt.csv");
    int count_non_zero = 0;
    for (std::size_t n1 = 0; n1 < d_secu.size1(); ++n1)
    {
        for (std::size_t n2 = 0; n2 < d_secu.size2(); ++n2)
        {
            if (d_opt(n1, n2) != 0) count_non_zero++;
            myfile << d_opt(n1, n2) << ",";
        }
        myfile << "\n";
    }
    myfile.close();

    cout << "#NON-ZERO D_OPT = " << count_non_zero / 2 << endl;
    */

    using domain_type = floe::domain::Domain;
    using time_mgr_type = floe::ope::TimeScaleManager<domain_type, decltype(detector)>;
    time_mgr_type time_mgr;
    domain_type domain;
    cout << "min DT : " << time_mgr.delta_t_secu(&domain, &detector);

    // SPECIFIC TEST FOR ID COUPLE
    size_t n1 = 122, n2 = 243;
    for (size_t i : {n1, n2})
        cout << i << " : " << "tau = " << detector.get_optim(i).tau() << ", center = " << detector.get_optim(i).global_disk().center << endl;
    
    // detect_step_3 analog
    auto const& opt1 = detector.get_optim(n1);
    auto const& opt2 = detector.get_optim(n2);

    // Intersection counter
    std::size_t cnt = 0;

    // Security and optimal distance
    value_type dist_s = std::numeric_limits<value_type>::max();
    value_type dist_o = dist_s;


    // Intersection finding
    for ( auto const& d1 : opt1.local_disks() )
    {
        for ( auto const& d2 : opt2.local_disks() )
        {
            const value_type dist = floe::collision::matlab::distance_circle_circle( d1, d2 );
            dist_o = std::min( dist_o, dist + d1.radius + d2.radius ); // unused

            if ( dist < 0 )
            {
                ++cnt;
            } else {
                dist_s = std::min( dist_s, dist + opt1.cdist() + opt2.cdist() );
            }
        }
    }
    cout << "STEP3 RESULT : " << dist_s << endl;
    // END SPECIFIC TEST FOR ID COUPLE


    // Create MAT file
    mat_t *matfp;
    matfp = Mat_CreateVer("io/test_dsecu.mat",NULL,MAT_FT_MAT5);
    if ( NULL == matfp ) {
    fprintf(stderr,"Error creating MAT file \"matfile5.mat\"!\n");
    cout << "EXIT_FAILURE" << endl;
    }

    // Write d_secu
    size_t N1 = d_secu.size1(), N2 = d_secu.size2();
    double data[N2][N1];
    for (std::size_t n1 = 0; n1 < N1; ++n1)
    {
        for (std::size_t n2 = 0; n2 < N1; ++n2)
        {
            data[n2][n1] = d_secu(n1, n2);
            // To place ghosts in same order than matlab for comparaison
            for (size_t k = 0; k<4; ++k)
                data[(k+1)*N1 + n2][n1] = d_secu(n1, N1 + 4 * n2 + k);
        }
    }

    matvar_t *matvar;
    size_t dims[2] = {N1,N2};
    matvar = Mat_VarCreate("dsecu",MAT_C_DOUBLE,MAT_T_DOUBLE,2,dims,data,0);
    if ( NULL == matvar ) {
    fprintf(stderr,"Error creating variable for ’dsecu’\n");
    } else {
    Mat_VarWrite(matfp,matvar,MAT_COMPRESSION_NONE);
    Mat_VarFree(matvar);
    }

    // Write d_opt
    // double dopt[N2][N1];
    for (std::size_t n1 = 0; n1 < N1; ++n1)
    {
        for (std::size_t n2 = 0; n2 < N1; ++n2)
        {
            data[n2][n1] = d_opt(n1, n2);
            // To place ghosts in same order than matlab for comparaison
            for (size_t k = 0; k<4; ++k)
                data[(k+1)*N1 + n2][n1] = d_opt(n1, N1 + 4 * n2 + k);
        }
    }

    matvar_t *matvar2;
    matvar2 = Mat_VarCreate("dopt",MAT_C_DOUBLE,MAT_T_DOUBLE,2,dims,data,0);
    if ( NULL == matvar2 ) {
    fprintf(stderr,"Error creating variable for ’dopt’\n");
    } else {
    Mat_VarWrite(matfp,matvar2,MAT_COMPRESSION_NONE);
    Mat_VarFree(matvar2);
    }

    // Close Matlab file
    Mat_Close(matfp);
    cout << "EXIT_SUCCESS" << endl;

    // cout optim's cdist
    for (size_t n = 0; n!= 10; ++n)
        cout << "mass : " << F.get_floes()[n].mass() << " area : " << F.get_floes()[n].area() << endl;

}
