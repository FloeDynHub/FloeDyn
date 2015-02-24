/*!
 * \file quadrature.hpp
 * \brief Quadrature base class
 * \author Roland Denis
 * \date November 2014
 */

#ifndef INTEGRATION_QUADRATURE_HPP_INCLUDED
#define INTEGRATION_QUADRATURE_HPP_INCLUDED

#include <type_traits>

namespace floe { namespace integration {
    
    /*! Base class for any quadrature method.
     *
     * It is a kind of abstract class using static polymorphism.
     *
     * \tparam M Quadrature method (for CRTP)
     * \tparam T Space fondamental type
     * \tparam D Space dimension
     */
    template <
        class M,        // Quadrature method (CRTP)
        typename T,     // Space fondamental type
        unsigned int D  // Space dimension
    >
    class Quadrature {
        public :
        
        // Type traits
        typedef T space_type;   //!< Space fondamental type
        enum { space_dim = D }; //!< Space dimension

        /*! Integration of a function on the reference triangle
         * \tparam Function Type of the function to integrate
         * \param f The function to integrate
         */
        template < typename Function >
        constexpr
        auto apply( const Function & f ) const 
            -> typename std::result_of<Function&(T,T)>::type
        {
            return static_cast<const M*>(this)->apply_impl(f);
        }
                
    };

}} // namespace floe::integration

#endif
