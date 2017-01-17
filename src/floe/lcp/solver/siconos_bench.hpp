/*!
 * \file lcp/solver/SiconosBench.hpp
 * \brief LCP solver
 * \author Quentin Jouet
 */

#ifndef LCP_SICONOS_BENCH_HPP
#define LCP_SICONOS_BENCH_HPP

#include <vector>
#include <array>
#include <string>
#include <siconos/LCP_Solvers.h>
#include <siconos/lcp_cst.h>
#include <chrono>
#include "floe/lcp/lcp.hpp"
#include "floe/lcp/solver/newton_min.hpp"


namespace floe { namespace lcp { namespace solver
{


struct SolverStat
{
    using lcp_type = floe::lcp::LCP<double>;
    SolverStat(LCP_SOLVER s_id, std::string s_name) : solver_id{s_id}, solver_name{s_name}, is_siconos{true}{
            global_options.verboseMode = 0;
        }
    SolverStat(int s_id, std::string s_name) : solver_id{s_id}, solver_name{s_name}, is_siconos{false}{
            global_options.verboseMode = 0;
        }

    int solver_id;
    std::string solver_name;
    SolverOptions base_options;
    NumericsOptions global_options;
    bool is_siconos;

    int nb_success{0};
    int nb_lcp{0};
    double compute_time{0.0};
    double error_on_success{0};

    bool solve(lcp_type& lcp){
        nb_lcp++;
        bool success{0};
        auto t_start = std::chrono::high_resolution_clock::now();
        if (is_siconos)
        {
            success = siconos_solve(lcp);
        } else {
            if (solver_id == 0){
                // auto lcp_copy = lcp; // newton_min seems to modify lcp (TODO test again)
                // success = newton_min(lcp_copy);
                // lcp.z = lcp_copy.z;
                success = newton_min(lcp);
            }
        }
        auto t_end = std::chrono::high_resolution_clock::now();
        compute_time += std::chrono::duration<double, std::milli>(t_end-t_start).count();
        if (success) {    
            nb_success++;
            auto err = LCP_error(lcp);
            if (err == err) { error_on_success += err; }//if (solver_id == 0) std::cout << "NEWTON** " << err << " "; }
            // std::cout << "** " << solver_name << " : " << lcp.z << "\n";
        }
        return success;
    }
    bool siconos_solve(lcp_type& lcp){
        auto A = createNumericsMatrixFromData(0, lcp.dim, lcp.dim, lcp.A.data().begin()); // first arg = storage (0-dense/1-sparse)
        auto LCP = new LinearComplementarityProblem();
        LCP->M = A;
        LCP->size = lcp.dim;
        LCP->q = lcp.q.data().begin();
        double* Zs{lcp.z.data().begin()};
        double* Ws{lcp.w.data().begin()};
        linearComplementarity_setDefaultSolverOptions(LCP, &base_options, solver_id);
        int info = linearComplementarity_driver(LCP, Zs,Ws, &base_options, &global_options);
        delete LCP;
        delete A;
        return (info == 0);
    }
    double success_rate() { return (double)nb_success/nb_lcp; }
    double avg_compute_time() { return (double)compute_time/nb_lcp; }
    double avg_success_error() { return (double)error_on_success/nb_success; }
    void print_result() { std::cout << solver_name << " : "
        << nb_success << " / " << nb_lcp << " (" << 100 * success_rate() << " %)"
        << " avg_time : " << avg_compute_time() << " ms"
        << ", avg_success_error : " << avg_success_error() << "\n";
    }

};


struct SiconosBench
{
    using lcp_type = floe::lcp::LCP<double>;
    SiconosBench() : solver_stat_list{
        {SICONOS_LCP_LEMKE, SICONOS_LCP_LEMKE_STR},
        // {SICONOS_LCP_NSQP, SICONOS_LCP_NSQP_STR}, // prints "***ql: matrix g was enlarged  1-times by unitmatrix"
        // {SICONOS_LCP_NEWTONMIN, SICONOS_LCP_NEWTONMIN_STR}, // nan solutions !!!
        {0, SICONOS_LCP_NEWTONMIN_STR}, // custom version
        // {SICONOS_LCP_NEWTON_FBLSA, SICONOS_LCP_NEWTON_FBLSA_STR}, // nan solutions !!!
        {SICONOS_LCP_RPGS, SICONOS_LCP_RPGS_STR},
        {SICONOS_LCP_AVI_CAOFERRIS, SICONOS_LCP_AVI_CAOFERRIS_STR},
        // {SICONOS_LCP_PIVOT, SICONOS_LCP_PIVOT_STR}, // prints "do_pivot_driftless :: pivot value too small 2.047463e-19"
    }, nb_lcp{0}, nb_global_success{0} {}

    ~SiconosBench() { print_results(); }

    std::vector<SolverStat> solver_stat_list;
    int nb_lcp;
    int nb_global_success;

    void run(lcp_type lcp){
        nb_lcp++;
        bool global_success = false;
        for (auto& stat : solver_stat_list)
        {
            if (stat.solve(lcp)) global_success = true; 
        }
        if (global_success) nb_global_success++;
    }
    double success_rate() { return (double)nb_global_success/nb_lcp; }
    void print_results(){
        std::cout << "*****SICONOS BENCH RESULTS*****\n";
        for (auto& stat : solver_stat_list)
        {
            stat.print_result();
        }
        std::cout << "GLOBAL : " << nb_global_success << " / " << nb_lcp << " (" << 100 * success_rate() << " %)\n";
    }
};

    // std::cout << "SICONOS INFO : " << info;

/* COMPLETE SOLVER LIST 
{SICONOS_LCP_LEMKE, SICONOS_LCP_LEMKE_STR},
        // {SICONOS_LCP_NSGS_SBM, SICONOS_LCP_NSGS_SBM_STR},
        // {SICONOS_LCP_PGS, SICONOS_LCP_PGS_STR},
        {SICONOS_LCP_CPG, SICONOS_LCP_CPG_STR}, // bad perf
        {SICONOS_LCP_LATIN, SICONOS_LCP_LATIN_STR}, // not adapted
        {SICONOS_LCP_LATIN_W, SICONOS_LCP_LATIN_W_STR}, // not adapted
        // {SICONOS_LCP_QP, SICONOS_LCP_QP_STR}, // not adapted
        {SICONOS_LCP_NSQP, SICONOS_LCP_NSQP_STR},
        {SICONOS_LCP_NEWTONMIN, SICONOS_LCP_NEWTONMIN_STR},
        // {SICONOS_LCP_NEWTON_FBLSA, SICONOS_LCP_NEWTON_FBLSA_STR},
        {SICONOS_LCP_PSOR, SICONOS_LCP_PSOR_STR}, //" Warning negative diagonal term  The local problem cannot be solved \n"
        {SICONOS_LCP_RPGS, SICONOS_LCP_RPGS_STR},
        {SICONOS_LCP_PATH, SICONOS_LCP_PATH_STR},  // not adapted (0%)
        {SICONOS_LCP_ENUM, SICONOS_LCP_ENUM_STR}, // bad perf
        {SICONOS_LCP_AVI_CAOFERRIS, SICONOS_LCP_AVI_CAOFERRIS_STR},
        {SICONOS_LCP_PIVOT, SICONOS_LCP_PIVOT_STR},
        // {SICONOS_LCP_BARD, SICONOS_LCP_BARD_STR},
        // {SICONOS_LCP_MURTY, SICONOS_LCP_MURTY_STR},
        // {SICONOS_LCP_NEWTON_MINFBLSA, SICONOS_LCP_NEWTON_MINFBLSA_STR}, // bad perf
        {SICONOS_LCP_PATHSEARCH, SICONOS_LCP_PATHSEARCH_STR}, // For testing only
        {SICONOS_LCP_PIVOT_LUMOD, SICONOS_LCP_PIVOT_LUMOD_STR}, // bad perf
        {SICONOS_LCP_GAMS, SICONOS_LCP_GAMS_STR}, // gams was not enabled at compile time!
        */

}}} // namespace floe::lcp::solver


#endif // LCP_SICONOS_BENCH_HPP