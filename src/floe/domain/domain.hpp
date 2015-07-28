/*!
 * \file domain/domain.hpp
 * \brief domain
 * \author Quentin Jouet
 */

#ifndef DOMAIN_DOMAIN_HPP
#define DOMAIN_DOMAIN_HPP
#include <iostream> // debug



namespace floe { namespace domain
{

/*! Domain
 *
 * Time management
 *
 */


class Domain
{

public:

    using real = VALUE_TYPE;

    // Default constructor
    Domain() : m_t{0}, m_delta_t{1}, m_delta_t_default{DT_DEFAULT}, m_last_out{-1} {}

    inline const real& time() const { return m_t; }
    inline void set_time(real t) { m_t = t; }
    inline const real& time_step() const {return m_delta_t; }
    inline real default_time_step() const { return m_delta_t_default; }
    inline void set_time_step(real delta_t ) { m_delta_t = delta_t; std::cout << "delta_t : " << delta_t << std::endl; }
    inline void update_time() { m_t += m_delta_t; }
    inline const real& last_out() const { return m_last_out; }
    inline void update_last_out() { m_last_out = m_t; }

private:

    real m_t;
    real m_delta_t;
    const real m_delta_t_default;

    real m_last_out;

};


}} // namespace floe::domain


#endif // DOMAIN_DOMAIN_HPP