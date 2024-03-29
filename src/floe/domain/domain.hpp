/*!
 * \file domain/domain.hpp
 * \brief Domain management
 * \author Quentin Jouet
 */

#ifndef DOMAIN_DOMAIN_HPP
#define DOMAIN_DOMAIN_HPP
#include <iostream> // debug
#include <cmath>



namespace floe { namespace domain
{

/*! Domain
 *
 * Time management
 *
 */

template<typename T>
class Domain
{

public:

    using real = T;

    //! Default constructor
    Domain() : m_t{0}, m_delta_t{1}, m_delta_t_default{10} {}

    inline const real& time() const { return m_t; }
    inline void set_time(real t) { m_t = t; }
    inline const real& time_step() const {return m_delta_t; }
    inline void set_time_step(real delta_t ) { m_delta_t = delta_t; }
    inline real default_time_step() const { return m_delta_t_default; }
    inline void set_default_time_step(real delta_t ) { m_delta_t_default = delta_t; }
    inline void update_time() { m_t += m_delta_t; }
    inline void rewind_time() { m_t -= m_delta_t; }

private:
    // time
    real m_t; //!< Simulation time
    real m_delta_t; //!< Time step
    real m_delta_t_default; //!< Default time step

};


}} // namespace floe::domain


#endif // DOMAIN_DOMAIN_HPP