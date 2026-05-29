/*!
 * \file lcp/LCP_manager.hpp
 * \brief LCP manager
 * \author Quentin Jouet
 */

#ifndef OPE_LCP_MANAGER_HPP
#define OPE_LCP_MANAGER_HPP

#include "floe/domain/time_scale_manager.hpp"
#include "../product/config/config_base.hpp" // types

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/blas.hpp>
#include <iostream> // debug
#include <array>
#include <algorithm>
#include <chrono> // test perf
#include <vector>
#include <unordered_map> // for m_subgraph_hash_count
#include <boost/graph/adjacency_list.hpp>

// saving matrix when lcp solver failed for further analysing
#include "H5Cpp.h"
#include "floe/utils/random.hpp"                    // for saving lcp statistic matrix 

// #include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/graph_utility.hpp>
// #include <boost/graph/filtered_graph.hpp>
// #include <boost/graph/connected_components.hpp>
// #include <boost/graph/copy.hpp>

// #include "floe/collision/contact_point.hpp"
// #include "floe/collision/floe_contact.hpp"
// #include "floe/collision/floe_vertex.hpp"


#ifdef _OPENMP
#include <omp.h>
#endif

namespace floe { namespace lcp
{

using namespace types;

/*! LCPManager
 *
 * Operator for Collision processing
 *
 */


template<typename TSolver>
class LCPManager
{

public:
    using solver_type = TSolver;
    using value_vector = boost::numeric::ublas::vector<real_type>;

    //! Constructor
    LCPManager(real_type epsilon) : m_solver{epsilon}, m_nb_lcp{0}, m_nb_lcp_success{0} {}

    //! Destructor
    ~LCPManager(){ 
        if (m_nb_lcp)
            std::cout   << "#TOTAL LCP solve: " << m_nb_lcp_success << "/" << m_nb_lcp << "(" << success_ratio() << "%) \n"
                        << "LCP_failed compression phase: " << m_nb_lcp_failed_stats[0] << 
                        ", LCP_failed decompression phase: " << m_nb_lcp_failed_stats[1] << "\n" <<
                        ", LCP_solved with solution maintaining the kinetic energy: " << m_nb_lcp_failed_stats[2] << "\n";
    }

    //! LCP solver accessor
    inline solver_type& get_solver() { return m_solver; }

    //! Solve collision represented by a contact graph
    template<typename TContactGraph>
    int solve_contacts(TContactGraph& contact_graph, real_type time, real_type dt);
    //! Get solving success ratio in percent
    double success_ratio(){ return (m_nb_lcp == 0)? 100 : 100 * (double)m_nb_lcp_success/m_nb_lcp; }

private:

    solver_type m_solver; //!< LCP Solver
    long        m_nb_lcp; //!< Total number of LCP managed
    long        m_nb_lcp_success; //!< Total number of LCP solving success
    long        m_nb_lcp_failed_stats[3]={0,0,0}; // LCP failed statistics: [nb LCP failed during compression phase,
    // nb LCP failed during decompression phase, nb LCP solved maintaining the kinetic energy in decompression phase].
    // Other LCP statistics are found in the matlab routine (see folder: io/outputs).
    
    double chrono_active_subgraph{0.0}; // test perf
    double max_chrono_active_subgraph{0.0}; // test perf
    
    std::unordered_map<std::size_t, int> m_subgraph_hash_count; //!< Map to track subgraph hash occurrences
    std::unordered_map<std::size_t, std::vector<typename types::point_type>> m_subgraph_hash_cp_impulses;
    std::unordered_map<std::size_t, std::vector<real_type>> m_subgraph_hash_floe_impulses;
    std::unordered_map<std::size_t, real_type> m_subgraph_map_timestep; //!< Map to track subgraph hash to time interval corresponding to the impulses stored

    //! Update floes state with LCP solution
    template<typename TContactGraph>
    void update_floes_state(TContactGraph& graph, const value_vector Sol, real_type time);
    //! Update floes impulses with LCP solution
    template<typename TContactGraph>
    void update_floes_impulses(TContactGraph& graph, real_type time);

    /*! \fn bool saving_contact_graph_in_hdf5(int lCP_count, std::size_t loop_count, std::size_t size_a_sub_graph, bool all_solved )
        \brief Saves information on the contact graph in the same file as LCP statistics.

        \param lcp_count the total number of treated LCP at the end of the "while loop"
        \param loop_count the total number of loop performed at the end of the "while loop"
        \param size_asubgraph the number of LCP generated from the first active sub-graphic
        \param all_soved a boolean to indicate whether all treated LCP during the "while loop" are solved

        \details the data are saved in the <a>nx6</a> matrix form, with <a>n</a> the number of sub-graph dealt with during the simulation.
        The output is a boolean to ensure data do not exceed a defined maximum storage.
        Each line of the matrix contains: the numerous of the last unsolved LCP stored in the h5 file (see LCPsolver::saving_LCP_in_hdf5),
        the numerous of the last solved LCP stored in the h5 file, lcp_count, loop_count, size_asubgraph and all_solved.  

        \warning Do not be confused with contact graph, sub-graph and active sub-graph (see LCPManager::solve_contacts). 
    */
    bool saving_contact_graph_in_hdf5(int lcp_count, std::size_t loop_count, std::size_t size_asubgraph, 
        bool all_solved, int contact_loop_stats[] );
};


template<typename T>
template<typename TContactGraph>
int LCPManager<T>::solve_contacts(TContactGraph& contact_graph, real_type time, real_type dt)
{

    auto const subgraphs = collision_subgraphs( contact_graph );
    int LCP_count=0, nb_success=0;
    int nb_lcp_failed_stats[3]={0,0,0};
    int optim_jam_start = 10;

    // Variable for OPTIMJAM guessed impulse
    typename types::point_type guessed_contact_impulse;

    // Map to track seen hashes in this call OPTIMJAM
    std::unordered_set<std::size_t> seen_hashes;

    const std::size_t limit_sup_loop_cnt    = 40000;//5000; // from Quentin: 1000 // matthias : 800
    const std::size_t limit_sup_nb_contact  =  800;//500; // from Quentin:   50

    // variables for contact informations:
    #ifdef LCPSTATS
        static bool end_recording = false;
    #endif

    for ( auto& subgraph : subgraphs )
    {
        // OPTIMJAM
        // is there active contact in the subgraph ?
        bool has_active_contact = false;
        for (auto e : boost::make_iterator_range(edges(subgraph))) {
            for (std::size_t i = 0; i < subgraph[e].size(); ++i) {
                if (subgraph[e][i].is_active()) {
                    has_active_contact = true;
                    break;
                }
            }
            if (has_active_contact) break;
        }
        if (!has_active_contact) continue; // Ne considérer que les subgraphs avec contacts actifs

        // Collect floe pointers to build a unique hash for the subgraph
        std::vector<std::uintptr_t> floe_ptrs;
        for (auto v : boost::make_iterator_range(vertices(subgraph))) {
            floe_ptrs.push_back(reinterpret_cast<std::uintptr_t>(subgraph[v].floe));
        }
        std::sort(floe_ptrs.begin(), floe_ptrs.end());

        // Hash: combine floe addresses and number of contacts
        /* V1
        std::size_t graph_hash = 0;
        for (auto ptr : floe_ptrs) {
            graph_hash ^= ptr + 0x9e3779b9 + (graph_hash << 6) + (graph_hash >> 2);
        }
        graph_hash ^= num_contacts(subgraph) + 0x9e3779b9 + (graph_hash << 6) + (graph_hash >> 2);
        */
        // V2 : use floe_ptr AND contact points (round position x and y to unit (1m))
        std::size_t graph_hash = 0;
        for (auto v : boost::make_iterator_range(vertices(subgraph))) {
            auto floe = subgraph[v].floe;
            std::uintptr_t ptr = reinterpret_cast<std::uintptr_t>(floe);
            graph_hash ^= ptr + 0x9e3779b9 + (graph_hash << 6) + (graph_hash >> 2);
        }
        for (auto e : boost::make_iterator_range(edges(subgraph))) {
            for (std::size_t i = 0; i < subgraph[e].size(); ++i) {
                auto contact = subgraph[e][i];
                auto pos = contact.frame.center();
                std::size_t contact_hash = 0;
                // TODO adapt resolution to the scale of the problem.
                double resolution = 100.0; // 1 unit in hash corresponds to 1m in position
                contact_hash ^= static_cast<std::size_t>(static_cast<std::int64_t>(std::llround(pos.x / resolution))) + 0x9e3779b9 + (contact_hash << 6) + (contact_hash >> 2);
                contact_hash ^= static_cast<std::size_t>(static_cast<std::int64_t>(std::llround(pos.y / resolution))) + 0x9e3779b9 + (contact_hash << 6) + (contact_hash >> 2);
                graph_hash ^= contact_hash + 0x9e3779b9 + (graph_hash << 6) + (graph_hash >> 2);
            }
        }


        seen_hashes.insert(graph_hash);

        // Update map: increment if exists, else set to 1
        auto it = m_subgraph_hash_count.find(graph_hash);
        if (it != m_subgraph_hash_count.end()) {
            it->second += 3;
        } else {
            m_subgraph_hash_count[graph_hash] = 1;
            m_subgraph_map_timestep[graph_hash] = 0;
            m_subgraph_hash_cp_impulses[graph_hash] = std::vector<typename types::point_type>(
                num_contacts(subgraph), typename types::point_type{0, 0}
            );
            m_subgraph_hash_floe_impulses[graph_hash] = std::vector<real_type>(
                num_vertices(subgraph), 0
            );
        }

        std::cout << "Subgraph hash: " << graph_hash
              << " | Floe count: " << floe_ptrs.size()
              << " | Contacts: " << num_contacts(subgraph)
              << " | Seen count: " << m_subgraph_hash_count[graph_hash]
              << std::endl;
        // END OPTIMJAM

        // OPTIMJAM : solve subgraph reapplaying impulses from previous timesteps

        if (m_subgraph_hash_count[graph_hash] >= 3 * optim_jam_start and m_subgraph_map_timestep[graph_hash] > dt * 10) {
            // int i = 0;
            // idea #2 compute floes average speed in the subgraph (excluding obstacles) and subtract it to each floe
            typename types::point_type total_speed = {0, 0};
            real_type total_mass = 0;
            real_type max_v2 = 0;
            for (auto v : boost::make_iterator_range(vertices(subgraph))) {
                if (!subgraph[v].floe->is_obstacle()) {
                    total_speed += subgraph[v].floe->state().speed * subgraph[v].floe->mass();
                    total_mass += subgraph[v].floe->mass();
                    auto v2 = subgraph[v].floe->kinetic_energy() / subgraph[v].floe->mass();
                    if (v2 > max_v2) {
                        max_v2 = v2;
                    }
                }
            }
            std::cout << "max v2 in subgraph: " << max_v2 << std::endl;
            if (total_mass > 0 && max_v2 < 0.01) { // TODO find good threshold for max_v2
                std::cout << "OPTIM JAM !" << std::endl;
                // if the subgraph is almost jammed, we set speed to 0 to avoid numerical issues in LCP solver
                // and speed up collision computing
                // auto avg_speed = total_speed / total_mass;
                // std::cout << "Average speed in subgraph: (" << avg_speed[0] << ", " << avg_speed[1] << ")" << std::endl;
                int k = 0;
                for (auto v : boost::make_iterator_range(vertices(subgraph))) {
                    if (!subgraph[v].floe->is_obstacle()) {
                        // subgraph[v].floe->state().speed -= 2 * avg_speed;
                        subgraph[v].floe->state().speed = {0, 0};
                        subgraph[v].floe->state().rot = 0;
                        subgraph[v].floe->state().set_jammed(true);
                    }
                    double guessed_impulse_norm = 1 * dt * m_subgraph_hash_floe_impulses[graph_hash][k] / m_subgraph_map_timestep[graph_hash];
                    subgraph[v].add_impulse_received(guessed_impulse_norm);
                    ++k;
                }
                // for (auto e : boost::make_iterator_range(edges(subgraph))) {
                //     for (std::size_t j = 0; j < subgraph[e].size(); ++j) {
                //         guessed_contact_impulse = 1 * dt * m_subgraph_hash_cp_impulses[graph_hash][i] / m_subgraph_map_timestep[graph_hash];
                //         // Print floe speeds before impulse
                //         // auto& floe_src = *subgraph[source(e, subgraph)].floe;
                //         // auto& floe_tgt = *subgraph[target(e, subgraph)].floe;
                //         // // idea #1 re-apply guessed impulse
                //         // subgraph[e][j].add_impulse_received(guessed_contact_impulse);
                //         // // get floe1 and floe 2 from edge and add impulse to them
                //         // subgraph[source(e, subgraph)].floe->apply_impulse(-subgraph[e][j].impulse_abs_frame(), subgraph[e][j].frame.center());
                //         // subgraph[target(e, subgraph)].floe->apply_impulse(subgraph[e][j].impulse_abs_frame(), subgraph[e][j].frame.center());
                //         // Set received impulses to floes (visu only)
                //         real_type impulse_norm = norm2(guessed_contact_impulse);
                //         subgraph[source(e, subgraph)].add_impulse_received(impulse_norm);
                //         subgraph[target(e, subgraph)].add_impulse_received(impulse_norm);
                //         *subgraph[e][j].floe_states_changed = true;
                //     }
                //     ++i;
                //     subgraph[e].mark_changed();
                // }
                mark_changed(subgraph); // indicates which floes have been modified
            }
        }

        // Active subgraph LCP strategy
        auto asubgraphs = active_subgraphs( subgraph );
        if (asubgraphs.size() == 0) {
            std::cout << "OptimJam success ! No active contact in the subgraph." << std::endl;
        }
        std::size_t loop_cnt    = 0;
        int loop_nb_success     = -1;
        bool active_quad_cut    = 0;

        // variables for contact informations:
        #ifdef LCPSTATS
            std::size_t size_a_sub_graph = asubgraphs.size();
        #endif
        bool all_solved = true;

        int contact_loop_stats[2]={0,0};    // number of contact points, indicator for be out of loop due to all success (1) or no success (0) 
                                            // or no enough iteration (2)
        contact_loop_stats[0] = static_cast<int>(num_contacts(subgraph));
        contact_loop_stats[1] = 1;

        // UTEST Save floes speed and rot before LCP
        std::vector<typename types::point_type> utest_before_speeds;
        std::vector<real_type> utest_before_rots;
        for (auto const v : boost::make_iterator_range(vertices(subgraph))) {
            utest_before_speeds.push_back(subgraph[v].floe->state().speed);
            utest_before_rots.push_back(subgraph[v].floe->state().rot);
        }

        while (asubgraphs.size() != 0
               && loop_cnt < std::min( 100 * num_contacts(subgraph), limit_sup_loop_cnt) // 60 * num_contacts(subgraph)
               && loop_nb_success !=0 )
        {
            loop_nb_success = 0;    // if no succes after one total path of contact graph, no use (nothing change)
                                    // to browse again the while loop. Thus we add "loop_nb_succes !=0".

            LCP_count += asubgraphs.size();

            for ( auto const& graph : asubgraphs ) // loop over the total number of active contact group
            {
                bool success;
                if (num_contacts(graph) > limit_sup_nb_contact){
                    std::cout << "Q4, nb contact:" << num_contacts(graph) << " \n";
                    auto qdct = quad_cut( graph );
                    std::cout << "nb quad_cut: " << qdct.size() << " \n";
                    active_quad_cut = 1;
                    for ( auto const& igraph : quad_cut( graph ) ){
                        auto Sol = m_solver.solve( igraph, success, nb_lcp_failed_stats );
                        mark_solved(igraph, success);
                        if (success) {++loop_nb_success;}
                        update_floes_state(igraph, Sol, time);
                    }
                } else {
                    auto Sol = m_solver.solve( graph, success, nb_lcp_failed_stats );
                    mark_solved(graph, success);
                    if (success) {++loop_nb_success;}
                    update_floes_state(graph, Sol, time); // updates the velocity of floes
                }

                mark_changed_parent(graph, subgraph); // indicates which floes have been modified
            }
            asubgraphs = active_subgraphs( subgraph ); // computes the new relative velocitoies from velocities of modified floes 
            nb_success += loop_nb_success;
            ++loop_cnt;

            if (loop_nb_success==0) {contact_loop_stats[1] = 0;}
        }
        if (asubgraphs.size() != 0)
        {
            all_solved = false;
            if (loop_nb_success!=0) {contact_loop_stats[1] = 2;}
            std::cout << "End of the while loop without resolution of all contacts!! nb contact: "<< num_contacts(subgraph) << "\n";
            for ( auto const& graph : asubgraphs ) mark_solved(graph, false);
        }
        if (loop_nb_success == 0) {
            std::cout << "No SUCCESS at all!! nb contact: " << num_contacts(subgraph) << " nb active subgraphs : " << asubgraphs.size() << "\n";
        }

        // Mat
        // Saving data on LCP:
        /*
         * Recovery of contact data (LCP_count, etc). Save in h5 file:
         */
        #ifdef LCPSTATS
            if (!end_recording && size_a_sub_graph!=0) {
                end_recording = saving_contact_graph_in_hdf5( LCP_count, loop_cnt, size_a_sub_graph, all_solved, contact_loop_stats );
            } 
        #endif
        // End saving data on LCP
        // EndMat
        // OPTIMJAM save impulses for this subgraph if solved and if optimization is active
        if (true) {
            // iter over edges and contacts to get impulses
            size_t i = 0;
            std::vector<std::pair<int, int>> contact_coords;
            std::vector<std::pair<std::uintptr_t, std::uintptr_t>> floe_ptr_pairs;
            std::vector<std::pair<double, double>> impulses_per_dt;
            for (auto e : boost::make_iterator_range(edges(subgraph))) {
            for (std::size_t j = 0; j < subgraph[e].size(); ++j) {
                auto impulse = subgraph[e][j].get_impulse_received();
                m_subgraph_hash_cp_impulses[graph_hash][i] += impulse;
                // Get contact point coordinates, round to int
                // auto pt = subgraph[e][j].frame.center();
                // contact_coords.emplace_back(static_cast<int>(std::round(pt[0])), static_cast<int>(std::round(pt[1])));
                // // Store floe pointer address pairs
                // auto src_ptr = reinterpret_cast<std::uintptr_t>(subgraph[source(e, subgraph)].floe);
                // auto tgt_ptr = reinterpret_cast<std::uintptr_t>(subgraph[target(e, subgraph)].floe);
                // floe_ptr_pairs.emplace_back(src_ptr, tgt_ptr);
                // // Store impulse/dt
                // impulses_per_dt.emplace_back(impulse[0] / dt, impulse[1] / dt);
                // ++i;
            }
            }
            int k = 0;
            for (auto v : boost::make_iterator_range(vertices(subgraph))) {
                auto impulse = subgraph[v].impulse();
                m_subgraph_hash_floe_impulses[graph_hash][k] += impulse;
                ++k;
            }
            
            // Print as python list
            // std::cout << "cp_coords.append([";
            // for (size_t k = 0; k < contact_coords.size(); ++k) {
            // std::cout << "(" << contact_coords[k].first << ", " << contact_coords[k].second << ")";
            // if (k + 1 != contact_coords.size()) std::cout << ", ";
            // }
            // std::cout << "])" << std::endl;

            // // Print floe pointer address pairs
            // std::cout << "floe_ptr_pairs.append([";
            // for (size_t k = 0; k < floe_ptr_pairs.size(); ++k) {
            // std::cout << "(" << floe_ptr_pairs[k].first << ", " << floe_ptr_pairs[k].second << ")";
            // if (k + 1 != floe_ptr_pairs.size()) std::cout << ", ";
            // }
            // std::cout << "])" << std::endl;

            // Print impulses/dt as python list
            // std::cout << "impulses_per_dt.append([";
            // for (size_t k = 0; k < impulses_per_dt.size(); ++k) {
            // std::cout << "(" << impulses_per_dt[k].first << ", " << impulses_per_dt[k].second << ")";
            // if (k + 1 != impulses_per_dt.size()) std::cout << ", ";
            // }
            // std::cout << "])" << std::endl;
            
            m_subgraph_map_timestep[graph_hash] += dt;
            // print m_subgraph_hash_cp_impulses for this hash
            // std::cout << "IMP " << m_subgraph_hash_floe_impulses[graph_hash][] << ", " << m_subgraph_map_timestep[graph_hash] << " = " << m_subgraph_hash_cp_impulses[graph_hash][0] / m_subgraph_map_timestep[graph_hash] << std::endl;
            // std::cout << "[";
            // for (auto const& imp : m_subgraph_hash_cp_impulses[graph_hash]) std::cout << "(" << imp[0] / m_subgraph_map_timestep[graph_hash] << ", " << imp[1] / m_subgraph_map_timestep[graph_hash] << "), ";
            // std::cout << "]," << std::endl;
            // std::cout << "[";
            // for (auto const& imp : m_subgraph_hash_floe_impulses[graph_hash]) std::cout << imp / m_subgraph_map_timestep[graph_hash] << ", ";
            // std::cout << "]," << std::endl;
            
        }
        // END OPTIMJAM

        
        // // UTEST save after speeds and rots
        // std::vector<typename types::point_type> utest_after_speeds;
        // std::vector<real_type> utest_after_rots;
        // for (auto const v : boost::make_iterator_range(vertices(subgraph))) {
        //     utest_after_speeds.push_back(subgraph[v].floe->state().speed);
        //     utest_after_rots.push_back(subgraph[v].floe->state().rot);
        // }
        
        // // UTEST now reset floe speeds and rots to before LCP and apply impulses stored in contact points to later check if we get the same result
        // size_t idx = 0;
        // for (auto v : boost::make_iterator_range(vertices(subgraph))) {
        //     subgraph[v].floe->state().speed = utest_before_speeds[idx];
        //     subgraph[v].floe->state().rot = utest_before_rots[idx];
        //     ++idx;
        // }

        // // UTEST Apply guessed impulses from previous steps to contact points and floes
        // size_t i = 0;
        // for (auto const& e : boost::make_iterator_range(edges(subgraph))) {
        //     for (std::size_t j = 0; j < subgraph[e].size(); ++j) {
        //         // Apply impulses to floe1 and floe2
        //         // if (!subgraph[source(e, subgraph)].floe->is_obstacle()) {
        //         //     std::cout << "Source floe area: " << subgraph[source(e, subgraph)].floe->area()
        //         //               << " -> received impulse (abs frame): x = " << -subgraph[e][j].impulse_abs_frame()[0]
        //         //               << ", y = " << -subgraph[e][j].impulse_abs_frame()[1]
        //         //               << " at contact point (" << subgraph[e][j].frame.center()[0]
        //         //               << ", " << subgraph[e][j].frame.center()[1] << ")" << std::endl;
        //         // }
        //         // if (!subgraph[target(e, subgraph)].floe->is_obstacle()) {
        //         //     std::cout << "Target floe area: " << subgraph[target(e, subgraph)].floe->area()
        //         //               << " -> received impulse (abs frame): x = " << subgraph[e][j].impulse_abs_frame()[0]
        //         //               << ", y = " << subgraph[e][j].impulse_abs_frame()[1]
        //         //               << " at contact point (" << subgraph[e][j].frame.center()[0]
        //         //               << ", " << subgraph[e][j].frame.center()[1] << ")" << std::endl;
        //         // }
        //         subgraph[source(e, subgraph)].floe->apply_impulse(-subgraph[e][j].impulse_abs_frame(), subgraph[e][j].frame.center());
        //         subgraph[target(e, subgraph)].floe->apply_impulse(subgraph[e][j].impulse_abs_frame(), subgraph[e][j].frame.center());
        //         ++i;
        //     }
        // }

        // // UTEST Check if we get the same result as after LCP
        // bool utest_success = true;
        // idx = 0;
        // for (auto v : boost::make_iterator_range(vertices(subgraph))) {
        //     auto speed_diff = subgraph[v].floe->state().speed - utest_after_speeds[idx];
        //     auto rot_diff = subgraph[v].floe->state().rot - utest_after_rots[idx];
        //     double speed_diff_norm = std::sqrt(speed_diff[0]*speed_diff[0] + speed_diff[1]*speed_diff[1]);
        //     if (speed_diff_norm > 1e-5 || std::abs(rot_diff) > 1e-5) {
        //     utest_success = false;
        //     std::cout << "UTEST FAILED for floe idx " << idx << ": speed diff = (" 
        //         << speed_diff[0] << ", " << speed_diff[1] << "), rot diff = " << rot_diff << std::endl;
        //     }
        //     ++idx;
        // }
        // if (utest_success) {
        //     std::cout << "UTEST SUCCESS: Floe states match after reapplying impulses." << std::endl;
        // }
        // End UTEST
        // Compute and print average distance between contact points in the subgraph
        double total_dist = 0.0;
        double min_dist = std::numeric_limits<double>::max();
        size_t contact_count = 0;
        for (auto const& e : boost::make_iterator_range(edges(subgraph))) {
            for (std::size_t j = 0; j < subgraph[e].size(); ++j) {
            double dist = subgraph[e][j].dist;
            total_dist += dist;
            if (dist < min_dist) min_dist = dist;
            ++contact_count;
            }
        }
        double sgec = 0.0;
        for (auto const v : boost::make_iterator_range(vertices(subgraph))) {
            sgec += subgraph[v].floe->kinetic_energy();
        }
        if (contact_count > 0) {
            double avg_dist = total_dist / contact_count;
            std::cout << "Average contact point distance in subgraph: " << avg_dist
                  << ", min: " << min_dist << " ; kinetic energy : " << sgec << std::endl;
        }
        
    }

    update_floes_impulses(contact_graph, time);
    m_nb_lcp += LCP_count;
    m_nb_lcp_success += nb_success;
    for (int i=0;i<3;++i){
        m_nb_lcp_failed_stats[i] += nb_lcp_failed_stats[i];
    }

    // Remove hashes not seen in this call OPTIMJAM
    for (auto it = m_subgraph_hash_count.begin(); it != m_subgraph_hash_count.end(); ) {
        if (seen_hashes.find(it->first) == seen_hashes.end() || it->second > 500) {
            std::size_t key = it->first;
            it->second -= 1;
            if (it->second <= 0 || it->second >= 500) {
                m_subgraph_hash_cp_impulses.erase(key);
                m_subgraph_hash_floe_impulses.erase(key);
                m_subgraph_map_timestep.erase(key);
                it = m_subgraph_hash_count.erase(it);
            } else {
                ++it;
            }
        } else {
            ++it;
        }
    }

    #ifndef MPIRUN
    if (LCP_count)
        std::cout << " #LCP solve: "<< nb_success << " / " << LCP_count << std::endl;
    #endif
    return nb_success;
}


template<typename T>
template<typename TContactGraph>
void LCPManager<T>::update_floes_state(TContactGraph& graph, const value_vector Sol, real_type time){

    for ( auto const v : boost::make_iterator_range( vertices(graph) ) )
    {
        graph[v].floe->state().speed = {Sol(3*v), Sol(3*v + 1)}; // fv_test
        graph[v].floe->state().rot = Sol(3*v + 2); // fv_test
    }
}

//! Update floes impulses from contact graph
template<typename T>
template<typename TContactGraph>
void LCPManager<T>::update_floes_impulses(TContactGraph& graph, real_type time){
    for ( auto const v : boost::make_iterator_range( vertices(graph) ) )
    {
        graph[v].floe->add_impulse(graph[v].impulse()); // fv_test
    }
    for ( auto const& edge : make_iterator_range( edges( graph ) ) )
    {
        for ( std::size_t i = 0; i < graph[edge].size(); ++i ) // iter over contacts
        {
            // Add contact point impulses to corresponding floes
            graph[source(edge, graph)].floe->add_contact_impulse(
                graph[edge][i].frame.center(), -graph[edge][i].impulse_abs_frame(),
                time);
            graph[target(edge, graph)].floe->add_contact_impulse(
                graph[edge][i].frame.center(), graph[edge][i].impulse_abs_frame(),
                time);
        }
    }
}


template<typename T>
bool LCPManager<T>::saving_contact_graph_in_hdf5(int LCP_count, std::size_t loop_count, std::size_t size_a_sub_graph,
 bool all_solved, int contact_loop_stats[] ){

    using namespace H5;                                         // proper way inside a function (to prevent extension to entire code)

    const H5std_string FILE_NAME("io/outputs/LCP_stats.h5");
    const H5std_string GROUP_NAME_I( "solved" ); // root group
    const H5std_string GROUP_NAME_II( "unsolved" ); // root group
    const H5std_string Last_Memb( "Last LCP" ); // to prevent similar LCP
    const H5std_string Contact_Graph_Info( "Contact Graph Info" );
    const hsize_t Max_storage_line = 50000;

    const int dim_stats = 8; // number of saved data on the contact loop.

    /*
     * Try block to detect exceptions raised by any of the calls inside it
     */
    try{
        /*
         * Turn off the auto-printing when failure occurs so that we can
         * handle the errors appropriately
         */
        Exception::dontPrint();
        /*
         * Create or Open a file.
         */
        H5File* file;
        try {
            file = new H5File( FILE_NAME, H5F_ACC_RDWR );
        } catch (...) {
            return false;
        }

        Group* M_solved = new Group(file->openGroup(GROUP_NAME_I));
        Group* M_unsolved = new Group(file->openGroup(GROUP_NAME_II));

        int last_lcp_uns[1];
        DataSet* dataset_LM_uns = new DataSet(M_unsolved->openDataSet( Last_Memb ));
        dataset_LM_uns->read( last_lcp_uns, PredType::NATIVE_INT );

        int last_lcp[1];
        DataSet* dataset_LM = new DataSet(M_solved->openDataSet( Last_Memb ));
        dataset_LM->read( last_lcp, PredType::NATIVE_INT );

        delete dataset_LM_uns;
        delete dataset_LM;

        DataSet* CGI;
        hsize_t dim_curr_cgi[2]={0,0};
        if (last_lcp_uns[0]==0 && last_lcp[0]==0) {
            delete M_unsolved;
            delete M_solved;
            delete file;
            return false;
        }
        else {
            try {
                CGI = new DataSet(file->openDataSet(Contact_Graph_Info));

                DataSpace fpsace_ind = CGI->getSpace();
                fpsace_ind.getSimpleExtentDims( dim_curr_cgi, NULL); // retrieves the current dimensions

                if (dim_curr_cgi[0] > Max_storage_line){
                    std::cout << "the maximum storage (" << Max_storage_line << ") for contact graph information is reached.\n";
                    delete M_unsolved;
                    delete M_solved;
                    delete CGI;
                    delete file; 
                    return true; 
                }
 
                hsize_t ext_size[2] = { dim_curr_cgi[0]+1, dim_curr_cgi[1] };
                CGI->extend( ext_size ); // extension with one new line 

            }
            catch (...) {
                // Creation of dataset to store information on contact graph and the "while loop":
                hsize_t dim_cg[2] = {1, dim_stats};
                hsize_t maxdims_cg[2] = {H5S_UNLIMITED, dim_stats}; // unlimited dataspace
                DataSpace space_cg( 2, dim_cg, maxdims_cg );
                DSetCreatPropList prop_cg; // Modify dataset creation property to enable chunking
                hsize_t chunk_dims_cg[2] = {1, dim_stats}; // with extendible dataset we cannot use contiguous but chunked dataset
                prop_cg.setChunk(2, chunk_dims_cg);
                CGI = new DataSet(file->createDataSet( Contact_Graph_Info, PredType::NATIVE_INT, space_cg, prop_cg ));
            }
        }

        int contact_stat[dim_stats];
        contact_stat[0] = last_lcp_uns[0];
        contact_stat[1] = last_lcp[0];
        contact_stat[2] = LCP_count;
        contact_stat[3] = static_cast<int>(loop_count);
        contact_stat[4] = static_cast<int>(size_a_sub_graph);
        contact_stat[5] =  (all_solved == true)? 1:0;
        contact_stat[6] = contact_loop_stats[0]; // number of contact points
        contact_stat[7] = contact_loop_stats[1]; // indicator for be out of loop due to all success (1) or no success (0) 
        // or no enough iteration (2)

        DataSpace fspace_cgi = CGI->getSpace();
        hsize_t dim_cgi[2] = {1,dim_stats}; 
        hsize_t offset_cgi[2] = {dim_curr_cgi[0], 0};
        fspace_cgi.selectHyperslab( H5S_SELECT_SET, dim_cgi, offset_cgi); // selection of the hyperslab
        DataSpace mspace_cgi( 2, dim_cgi );
        CGI->write(contact_stat, PredType::NATIVE_INT, mspace_cgi, fspace_cgi); // write in the hyperslab
  
            
        delete M_unsolved;
        delete M_solved;
        delete CGI;
        delete file;
    }  // end of try block
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
    error.printErrorStack();
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
    error.printErrorStack();
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
    error.printErrorStack();
    }
    return false;
}


}} // namespace floe::lcp


#endif // OPE_LCP_MANAGER_HPP