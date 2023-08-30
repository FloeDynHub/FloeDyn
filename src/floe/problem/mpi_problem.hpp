/*!
 * \file floe/problem/problem.hpp
 * \brief Virtual Smooth mpi problem (parent of worker and master specializations)
 * \author Quentin Jouet
 */

#ifndef PROBLEM_MPI_PROBLEM_HPP
#define PROBLEM_MPI_PROBLEM_HPP

#include "floe/io/mpi_utils.hpp"

namespace floe { namespace problem
{

template <
    typename TProblem
>
class MPIProblem : public TProblem
{
public:
    using base_class = TProblem;
    using real_type = typename TProblem::real_type;
    using message_type = floe::io::InterProcessMessage<real_type>;
    using mpi_terminal_type = floe::io::MPIWrapper;

    //! Default constructor
    MPIProblem(real_type epsilon, int OBL_status, int buffer_size=1e6) :
        base_class(epsilon, OBL_status), m_mpi_term{buffer_size} {}

    void load_config(std::string const& filename) override {
        base_class::load_config(filename);
        // Store floes ids
        for (std::size_t i=0; i<this->get_floe_group().nb_floes(); ++i){
            this->get_floe_group().get_floes()[i].set_id(i);
        }
    }

    //! Solver of the problem (main method)
    virtual void solve(
        real_type end_time,
        real_type dt_default,
        real_type out_step = 0,
        bool reset = true,
        bool fracture = false,
        bool melting = false) override = 0;
    virtual mpi_terminal_type& mpi() { return m_mpi_term; }

private:
    //! MPI terminal
    mpi_terminal_type m_mpi_term;

};

}} // namespace floe::problem


#endif // PROBLEM_MPI_PROBLEM_HPP