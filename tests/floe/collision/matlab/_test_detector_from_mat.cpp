/*!
 * \file test_detector_from_mat.cpp
 * \brief Test collision detection from MATLAB list_so variable.
 * \author Roland Denis, Quentin Jouet
 *
 * The MATLAB file must be given as parameter of the program.
 */
#include "../tests/catch.hpp"
#include <iostream>

#include "floe/geometry/geometry.hpp"

#include "floe/io/matlab/list_so_to_floes.hpp"
#include "floe/io/matlab/list_so_import.hpp"
#include "floe/io/matlab/list_so.hpp"

#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"

#include "floe/collision/matlab/detector.h"

#include "floe/lcp/builder/graph_to_lcp.hpp"
// #include "floe/lcp/solver/lexicolemke_eigen.hpp"
#include "floe/lcp/solver/lexicolemke.hpp"
#include "floe/lcp/solver/lemke_eigen.hpp"
#include <boost/numeric/ublas/io.hpp>

TEST_CASE( "Test Detector from matlab", "[collision]" ) {
    using namespace std;
    using namespace floe::io::matlab;
    namespace fg = floe::geometry;
    namespace ff = floe::floes;
    
    using real = double;
    using TStaticFloe = ff::StaticFloe<real>;
    using TKinematicFloe = ff::KinematicFloe<TStaticFloe>;
    using TDetector = floe::collision::matlab::MatlabDetector<TKinematicFloe>;

    std::string mat_file_name = "tests/floe/io/matlab/r1day_set_up_250sm_sz_60_list_so_350_str.mat";

    MatlabListSolid<double> list_so;
    cout << "Reading \"" << mat_file_name << "\" ... " << endl;
    auto const& now = chrono::high_resolution_clock::now;
    auto t_start = now();
    read_list_so_from_file( mat_file_name, list_so);
    auto t_end = now();
    cout << "Chrono : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;

    cout << "Importing floes ... " << endl;
    std::vector<TKinematicFloe> floe_list;
    t_start = now();
    list_so_to_floes( list_so, floe_list );
    t_end = now();
    cout << "Chrono : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;

    // Initializing detector
    cout << "Linking floes to detector ..." << endl;
    t_start = now();
    TDetector detector;
    for (auto& floe : floe_list)
    {
        floe.update(); // TODO Avoid that ! (move doesn't move internal pointers)
        detector.push_back(&floe);
    }
    t_end = now();
    cout << "Chrono : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;

    // Update collisions
    cout << "Updating optimizers and detecting contacts ..." << endl;
    t_start = now();
    for (std::size_t i = 0; i < 1; ++i)
        detector.update();
    t_end = now();
    cout << "Chrono : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;

    cout << "Informations:" << endl;
    cout << "\t" << floe_list.size() << " floes." << endl;
    cout << "\t" << detector.num_points() << " points." << endl;
    cout << "\t" << detector.num_local_disks() << " local disks." << endl;

    const std::size_t num = num_contacts(detector.contact_graph());
    cout << "\t" << num << " contacts (" << (num*100./detector.num_points()) << "%)" << endl;
    cout << endl;

    /*
    cout << "Contact list:" << endl;
    //for ( auto const& contact : detector.contact_list() )
    for ( std::size_t i = 0; i < detector.contact_list().size(); ++i )
    {
        auto const& contact = detector.contact_list()[i];
        cout    << i+1 << ") "
                // << contact.n1+1 << ":" << contact.n2+1 << " ; "
                << "center = " << fg::dsv( contact.frame.center() ) << " ; "
                << "u = " << fg::dsv( contact.frame.u() ) << " ; "
                << "v = " << fg::dsv( contact.frame.v() ) << endl;
    }
    */

    cout << "Finding connected components in the contact graph ..." << endl;
    t_start = now();
    auto const subgraphs = collision_subgraphs( detector.contact_graph() );
    t_end = now();
    cout << "Chrono : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;
    cout << "\t" << subgraphs.size() << " connected components in the contact graph:" << endl << endl;

    for ( std::size_t i = 0 ; i < subgraphs.size(); ++i )
    {
        auto const& graph = subgraphs[i];
        cout << num_vertices(graph) << " floes in the component n°" << i << " (" << num_contacts(graph) << " contacts)" << endl;

        t_start = now();
        auto const asubgraphs = active_subgraphs( graph );
        t_end = now();
        cout << "Chrono : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;
        cout << "\twith " << asubgraphs.size() << " active sub-components ( ";

        size_t contact_cnt = 0;
        for ( auto const& g : asubgraphs )
        {
            contact_cnt += num_contacts(g);
            cout << num_vertices(g) << "(" << num_contacts(g) << ") ";
        }
        cout << ")" << endl;

        cout << "\t" << contact_cnt << " contacts taken into account." << endl << endl;
    }

    cout << "First floe informations:" << endl;
    cout << "\tdensity = " << floe_list[0].get_density() << endl;
    cout << "\tmass = " << floe_list[0].mass() << endl;
    cout << "\tarea = " << floe_list[0].area() << endl;
    cout << "\tdenom = " << floe_list[0].moment_cst() << endl;
    cout << endl;

    cout << "Preparation of one small LCP ..." << endl;
    {
    auto const test_graph = active_subgraphs( subgraphs[0] )[0];
    t_start = now();
    floe::lcp::builder::GraphLCP<real, decltype(test_graph)> graph_lcp( test_graph );
    t_end = now();
    cout << "Chrono : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;

    t_start = now();
    auto lcp = graph_lcp.getLCP();
    t_end = now();
    cout << "Chrono : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;
    cout << "\t of size " << lcp.A.size1() << "x" << lcp.A.size2() << endl;

    /*
    cout << "M = " << graph_lcp.M << endl;
    cout << "J = " << graph_lcp.J << endl;
    cout << "D = " << graph_lcp.D << endl;
    cout << "A = " << lcp.A << endl;
    cout << "q = " << lcp.q << endl;
    */

    cout << "Solving it ... " << flush;
    t_start = now();
    const bool success = floe::lcp::solver::lexicolemke(lcp);
    t_end = now();
    cout << success << endl;
    cout << "Chrono : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;
    cout << "\tErr = " << LCP_error(lcp) << endl;
    cout << "z = " << lcp.z << endl;
    
    cout << endl;
    }

    cout << "Preparation of one big LCP ..." << endl;
    {
    auto const test_graph = detector.contact_graph();
    t_start = now();
    floe::lcp::builder::GraphLCP<real, decltype(test_graph)> graph_lcp( test_graph );
    t_end = now();
    cout << "Chrono : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;

    t_start = now();
    auto lcp = graph_lcp.getLCP();
    t_end = now();
    cout << "Chrono : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;
    cout << "\t of size " << lcp.A.size1() << "x" << lcp.A.size2() << endl;
    cout << endl;
    }

    cout << "Finding all active connected components and building all LCPs ..." << endl;
    {
        size_t lcp_cnt = 0;
        size_t pt_cnt = 0;
        t_start = now();
        for ( size_t i = 0; i < subgraphs.size(); ++i )
        {
            auto const asubgraphs = active_subgraphs( subgraphs[i] );
            for ( auto const& graph : asubgraphs )
            {
                floe::lcp::builder::GraphLCP<real, decltype(graph)> graph_lcp( graph );
                auto lcp = graph_lcp.getLCP();
                ++lcp_cnt;
                pt_cnt += lcp.A.size1() * lcp.A.size2();
            }
        }
        t_end = now();
        cout << endl;
        cout << "Chrono : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;
        cout << "\t" << lcp_cnt << " LCPs with " << pt_cnt << " values." << endl;
        cout << endl;
    }

    cout << "Finding all active connected components, building all LCPs and solving it ..." << endl;
    {
        size_t lcp_cnt = 0;
        size_t pt_cnt = 0;
        t_start = now();
        for ( size_t i = 0; i < subgraphs.size(); ++i )
        {
            auto const asubgraphs = active_subgraphs( subgraphs[i] );
            for ( auto const& graph : asubgraphs )
            {
                floe::lcp::builder::GraphLCP<real, decltype(graph)> graph_lcp( graph );
                auto lcp = graph_lcp.getLCP();

                ++lcp_cnt;
                pt_cnt += lcp.A.size1() * lcp.A.size2();
                
                auto t_start2 = now();
                const bool success = floe::lcp::solver::lexicolemke(lcp);
                // const bool success = floe::lcp::solver::lemke(lcp);
                // const bool success = (floe::lcp::solver::lemke(lcp) || floe::lcp::solver::lexicolemke(lcp))
                auto t_end2 = now();
                cout << success << " | " ;//<< success2 << " | ";
                cout << "\t" << lcp.A.size1() << "x" << lcp.A.size2() << " : Err = " << LCP_error(lcp) << " ; "
                     << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms";
            }
        }
        t_end = now();
        cout << endl;
        cout << "Chrono : " << chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;
        cout << "\t" << lcp_cnt << " LCPs with " << pt_cnt << " values." << endl;
        cout << endl;
    }

}

