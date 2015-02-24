/*!
 * \file floe/floes/floe_exception.hpp
 * \brief Exception associated to floes.
 * \author Roland Denis
 */

#ifndef FLOE_FLOES_FLOE_EXCEPTION_HPP
#define FLOE_FLOES_FLOE_EXCEPTION_HPP

#include <exception>
#include <string>

namespace floe { namespace floes
{

//! Floe exception
class FloeException
    : std::exception
{
public:
    /*! Create an exception with message
     * \param message The message.
     */
    explicit FloeException( std::string const& message = "" ) : m_message(message) {}
    
    //! Return the message associated to this exception.
    virtual const char* what() const throw() { return m_message.c_str(); }

    //! Destructor.
    virtual ~FloeException() {}

private:

    std::string m_message; 
};


}} // namespace floe::floes

#endif // FLOE_FLOES_FLOE_EXCEPTION_HPP
