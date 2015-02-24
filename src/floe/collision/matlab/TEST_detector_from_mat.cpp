/*!
 * \file TEST_detector_from_mat.cpp
 * \brief Test collision detection from MATLAB list_so variable.
 * \author Roland Denis
 *
 * The MATLAB file must be given as parameter of the program.
 */

#include <iostream>

#include <boost/timer/timer.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "floe/geometry/geometry.hpp"

#include "floe/io/matlab/list_so_to_floes.hpp"
#include "floe/io/matlab/list_so_import.hpp"
#include "floe/io/matlab/list_so.hpp"

#include "floe/floes/static_floe.hpp"
#include "floe/floes/kinematic_floe.hpp"

#include "floe/collision/matlab/detector.hpp"

#include "floe/lcp/builder/graph_to_lcp.hpp"

int main( int argc, char* argv[] )
{
    using namespace std;
    using namespace floe::io::matlab;
    namespace fg = floe::geometry;
    namespace ff = floe::floes;
    
    using real = double;
    using TStaticFloe = ff::StaticFloe<real>;
    using TKinematicFloe = ff::KinematicFloe<TStaticFloe>;
    using TDetector = floe::collision::matlab::MatlabDetector<TKinematicFloe>;

    if ( argc < 2 )
    {
        cout << "Usage: " << argv[0] << " <matlab_file_name>" << endl;
        return 1;
    }

    MatlabListSolid<double> list_so;
    cout << "Reading \"" << argv[1] << "\" ... " << endl;
    boost::timer::cpu_timer timer;
    read_list_so_from_file( argv[1], list_so);
    timer.stop();
    cout << "\t" << timer.format() << endl;

    cout << "Importing floes ... " << endl;
    timer.start();
    auto floe_list = list_so_to_floes<TKinematicFloe>( list_so );
    timer.stop();
    cout << "\t" << timer.format() << endl;

    // Initializing detector
    cout << "Linking floes to detector ..." << endl;
    timer.start();
    TDetector detector;
    for (auto& floe_ptr : floe_list)
        detector.push_back(floe_ptr);
    timer.stop();
    cout << "\t" << timer.format() << endl;

    // Update collisions
    cout << "Updating optimizers and detecting contacts ..." << endl;
    timer.start();
    for (std::size_t i = 0; i < 1; ++i)
        detector.update();
    timer.stop();
    cout << "\t" << timer.format() << endl;

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
    timer.start();
    auto const subgraphs = collision_subgraphs( detector.contact_graph() );
    timer.stop();
    cout << "\t" << timer.format();
    cout << "\t" << subgraphs.size() << " connected components in the contact graph:" << endl << endl;

    for ( std::size_t i = 0 ; i < subgraphs.size(); ++i )
    {
        auto const& graph = subgraphs[i];
        cout << num_vertices(graph) << " floes in the component n°" << i << " (" << num_contacts(graph) << " contacts)" << endl;

        timer.start();
        auto const asubgraphs = active_subgraphs( graph );
        timer.stop();
        cout << "\t" << timer.format();
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
    cout << "\tdensity = " << floe_list[0]->get_density() << endl;
    cout << "\tmass = " << floe_list[0]->mass() << endl;
    cout << "\tarea = " << floe_list[0]->area() << endl;
    cout << "\tdenom = " << floe_list[0]->moment_cst() << endl;
    cout << endl;

    cout << "Preparation of one small LCP ..." << endl;
    {
    auto const test_graph = active_subgraphs( subgraphs[0] )[0];
    timer.start();
    floe::lcp::builder::GraphLCP<real, decltype(test_graph)> graph_lcp( test_graph );
    timer.stop();
    cout << "\t" << timer.format();

    timer.start();
    auto const A = graph_lcp.getLCP();
    timer.stop();
    cout << "\t" << timer.format();
    cout << "\t of size " << A.size1() << "x" << A.size2() << endl;
    
    cout << "M = " << graph_lcp.M << endl;
    cout << "J = " << graph_lcp.J << endl;
    cout << "D = " << graph_lcp.D << endl;
    cout << "A = " << A << endl;
    
    }

    cout << "Preparation of one big LCP ..." << endl;
    {
    auto const test_graph = detector.contact_graph();
    timer.start();
    floe::lcp::builder::GraphLCP<real, decltype(test_graph)> graph_lcp( test_graph );
    timer.stop();
    cout << "\t" << timer.format();

    timer.start();
    auto const A = graph_lcp.getLCP();
    timer.stop();
    cout << "\t" << timer.format();
    cout << "\t of size " << A.size1() << "x" << A.size2() << endl;
    }

    cout << "Finding all active connected components and building all LCPs ..." << endl;
    {
        size_t lcp_cnt = 0;
        size_t pt_cnt = 0;
        timer.start();
        for ( size_t i = 0; i < subgraphs.size(); ++i )
        {
            auto const asubgraphs = active_subgraphs( subgraphs[i] );
            for ( auto const& graph : asubgraphs )
            {
                floe::lcp::builder::GraphLCP<real, decltype(graph)> graph_lcp( graph );
                auto const A = graph_lcp.getLCP();
                ++lcp_cnt;
                pt_cnt += A.size1()*A.size2();

            }
        }
        timer.stop();
        cout << "\t" << timer.format();
        cout << "\t" << lcp_cnt << " LCPs with " << pt_cnt << " values." << endl;
    }

    


    // Freeing memory
    for ( auto& floe_ptr : floe_list )
        delete floe_ptr;
    

    return 0;
}

