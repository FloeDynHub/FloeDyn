/*!\file floe/arithmetic/container_operators.hpp
 * \brief Bunch of generic operators for containers.
 * \author Roland Denis
 *
 * Originates from FIRST project.
 *
 * \namespace container_operators
 * \brief Generic operators for containers.
 */

// #ifndef FLOE_ARITHMETIC_CONTAINER_OPERATORS_HPP
// #define FLOE_ARITHMETIC_CONTAINER_OPERATORS_HPP

// IDEE : ça devrait pouvoir être adaptable aux tuple !!!
// 		* faut déjà remplacer l'opérateur [] par get<id>
// 		* les opérations de réduction se sont plus compatibles (mais compilation OK si on ne les appelle pas)
//		* comment se comportent les lambdas avec un type dépendant de l'indice ? (template & lambda ?)

#include <cmath>
#include <cstddef>     // std::size_t
#include <type_traits> // std::enable_if
#include <ostream>

// namespace floe { namespace arithmetic {

namespace container_operators {

using std::size_t;

/*! Standard implementation to get the size, valueType and construct & fill new object.
 *
 * Work perfectly with std::array ;)
 */
template< typename T >
struct size_impl {
	static constexpr size_t size() {
		return T().size();
	}
};

template< typename T,
	class = typename std::enable_if< ! std::is_fundamental<T>::value >::type >
struct valueType
{
	typedef typename T::valueType type; //!< We can probably also use decltype(T()[0])
};

template< typename T >
struct cstr_impl {
	template< typename... V >
	static constexpr T from_values( V... values ) {
		return T({{ values... }});
	}
};

/*
 * Switch to activate the following operators for specific type
 */

struct yes { static const bool value = true;  };
struct no  { static const bool value = false; };

//! To be specialized with a type to activate
template< class T >
struct is_activated : no {}; // By default, desactivated for all types

template< class T, class Return >
struct activate : std::enable_if< is_activated<T>::value , Return > {};

/*
 * Some internal tools
 */
namespace tools {
/* 
 * Generate a templated integer sequence from 0 to N-1
 * Source : http://stackoverflow.com/questions/7858817/unpacking-a-tuple-to-call-a-matching-function-pointer
 */

template< size_t... >
struct seq { };

template< size_t N, size_t... S >
struct gen_seq : gen_seq< N-1, N-1, S... > { };

template< size_t... S >
struct gen_seq<0, S...> {
  typedef seq< S... > type;
};

/*
 * Apply a index-dependant function to a container
 */
template< class T, class Function, size_t... I >
constexpr T map_impl( Function fn, seq<I...> ) {
		return cstr_impl<T>::from_values( fn(I)... );
}

template< class T, class Function >
constexpr T map( Function fn ) {
        return map_impl<T>( fn, typename gen_seq< size_impl<T>::size() >::type() );
}

/*
 * Apply a compound operator to a container
 */

template < class... T >
inline
void compound_dummy( T... )
{}

template < class Function, size_t... I >
inline
void compound_impl( Function fn, seq<I...> )
{
    compound_dummy(fn(I)...) ;
}

template < class T, class Function >
inline
void compound( Function fn )
{
    compound_impl( fn, typename gen_seq< size_impl<T>::size() >::type() );
}

/*
 * Apply a reduction operator to a container
 * TODO: at compile-time
 */
template< class T, class Operator >
typename valueType<T>::type
reduce( const T & lhs, Operator op, typename valueType<T>::type init ) {
        for (size_t i = 0 ; i < size_impl<T>::size() ; ++i)
                init = op(init, lhs[i]);
        return init;
}

	
} // namespace tools

} // namespace container_operators

// }} // namespace floe::arithmetic

using namespace container_operators;

// Standard functions
using std::sqrt;
using std::abs;
using std::min;
using std::max;
using std::pow;
using std::cos;
using std::sin;
using std::tan;


/*
 * Binary operators
 */

// Addition
template< class V, class T >
constexpr
typename activate<T, T>::type
operator+( V lhs, const T & rhs ) {
	return tools::map<T>( [&] (size_t i) { return lhs + rhs[i]; } );
}

template< class T, class V >
constexpr
typename activate<T, T>::type
operator+( const T & lhs, V rhs ) {
	return rhs + lhs;	
}

template< class T >
constexpr
typename activate<T, T>::type
operator+( const T & lhs, const T & rhs ) {
	return tools::map<T>( [&] (size_t i) { return lhs[i] + rhs[i]; } );
}

// Subtraction
template< class V, class T >
constexpr
typename activate<T, T>::type
operator-( V lhs, const T & rhs ) {
	return tools::map<T>( [&] (size_t i) { return lhs - rhs[i]; } );
}

template< class T, class V >
constexpr
typename activate<T, T>::type
operator-( const T & lhs, V rhs ) {
	return tools::map<T>( [&] (size_t i) { return lhs[i] - rhs; } );	
}

template< class T >
constexpr
typename activate<T, T>::type
operator-( const T & lhs, const T & rhs ) {
	return tools::map<T>( [&] (size_t i) { return lhs[i] - rhs[i]; } );
}

// Multiplication
template< class V, class T >
constexpr
typename activate<T, T>::type
operator*( V lhs, const T & rhs ) {
	return tools::map<T>( [&] (size_t i) { return lhs * rhs[i]; } );
}

template< class T, class V >
constexpr
typename activate<T, T>::type
operator*( const T & lhs, V rhs ) {
	return tools::map<T>( [&] (size_t i) { return lhs[i] * rhs; } );	
}

template< class T >
constexpr
typename activate<T, T>::type
operator*( const T & lhs, const T & rhs ) {
	return tools::map<T>( [&] (size_t i) { return lhs[i] * rhs[i]; } );
}

// Division
template< class V, class T >
constexpr
typename activate<T, T>::type
operator/( V lhs, const T & rhs ) {
	return tools::map<T>( [&] (size_t i) { return lhs / rhs[i]; } );
}

template< class T, class V >
constexpr
typename activate<T, T>::type
operator/( const T & lhs, V rhs ) {
	return tools::map<T>( [&] (size_t i) { return lhs[i] / rhs; } );	
}

template< class T >
constexpr
typename activate<T, T>::type
operator/( const T & lhs, const T & rhs ) {
	return tools::map<T>( [&] (size_t i) { return lhs[i] / rhs[i]; } );
}

// Power
template< class V, class T >
constexpr
typename activate<T, T>::type
pow( V lhs, const T & rhs ) {
	return tools::map<T>( [&] (size_t i) { return pow( lhs, rhs[i] ); } );
}

template< class T, class V >
constexpr
typename activate<T, T>::type
pow( const T & lhs, V rhs ) {
	return tools::map<T>( [&] (size_t i) { return pow( lhs[i], rhs ); } );	
}

template< class T >
constexpr
typename activate<T, T>::type
pow( const T & lhs, const T & rhs ) {
	return tools::map<T>( [&] (size_t i) { return pow( lhs[i], rhs[i] ); } );
}

// Element-Wise min & max
template< class T >
constexpr
typename activate<T, T>::type
ew_min( const T & lhs, const T & rhs ) {
	return tools::map<T>( [&] (size_t i) { return min( lhs[i], rhs[i] ); } );
}

template< class T >
constexpr
typename activate<T, T>::type
ew_max( const T & lhs, const T & rhs ) {
	return tools::map<T>( [&] (size_t i) { return max( lhs[i], rhs[i] ); } );
}

// Compound addition
template< class T >
inline
typename activate<T, T&>::type
operator+=( T & lhs, const T & rhs )
{
    tools::compound<T>( [&] (size_t i) { return lhs[i] += rhs[i]; } );
    return lhs;
}

template< class T, class V >
inline
typename activate<T, T&>::type
operator+=( T & lhs, const V rhs )
{
    tools::compound<T>( [&] (size_t i) { return lhs[i] += rhs; } );
    return lhs;
}

// Compound substraction
template< class T >
inline
typename activate<T, T&>::type
operator-=( T & lhs, const T & rhs )
{
    tools::compound<T>( [&] (size_t i) { return lhs[i] -= rhs[i]; } );
    return lhs;
}

template< class T, class V >
inline
typename activate<T, T&>::type
operator-=( T & lhs, const V rhs )
{
    tools::compound<T>( [&] (size_t i) { return lhs[i] -= rhs; } );
    return lhs;
}

// Compound multiplication
template< class T >
inline
typename activate<T, T&>::type
operator*=( T & lhs, const T & rhs )
{
    tools::compound<T>( [&] (size_t i) { return lhs[i] *= rhs[i]; } );
    return lhs;
}

template< class T, class V >
inline
typename activate<T, T&>::type
operator*=( T & lhs, const V rhs )
{
    tools::compound<T>( [&] (size_t i) { return lhs[i] *= rhs; } );
    return lhs;
}

// Compound division
template< class T >
inline
typename activate<T, T&>::type
operator/=( T & lhs, const T & rhs )
{
    tools::compound<T>( [&] (size_t i) { return lhs[i] /= rhs[i]; } );
    return lhs;
}

template< class T, class V >
inline
typename activate<T, T&>::type
operator/=( T & lhs, const V rhs )
{
    tools::compound<T>( [&] (size_t i) { return lhs[i] /= rhs; } );
    return lhs;
}


/*
 * Unary operators
 */
#define FAST_UNARY( op_name , expression ) \
template< class T >                        \
constexpr                                  \
typename activate<T, T>::type              \
op_name ( const T & arg ) { return tools::map<T>( [&] (size_t i) { return expression ; } ); }

FAST_UNARY( operator+ ,  arg[i] )
FAST_UNARY( operator- , -arg[i] )

/*
 * Mathematical functions
 */
FAST_UNARY( cos  , cos(arg[i])  )
FAST_UNARY( sin  , sin(arg[i])  )
FAST_UNARY( tan  , tan(arg[i])  )
FAST_UNARY( sqrt , sqrt(arg[i]) )
FAST_UNARY( abs  , abs(arg[i])  )

/*
 * Reduction operations
 */
template< class T >
typename activate<T, typename valueType<T>::type >::type
sum( const T & arg ) {
	return tools::reduce( arg, std::plus<typename valueType<T>::type>(), typename valueType<T>::type(0) );
}

template< class T >
typename activate<T, typename valueType<T>::type >::type
prod( const T & arg ) {
	return tools::reduce( arg, std::multiplies<typename valueType<T>::type>(), typename valueType<T>::type(1) );
}

template< class T >
typename activate<T, typename valueType<T>::type >::type
norm2( const T & arg ) {
	return sqrt(sum(arg * arg));
}

template< class T >
typename activate<T, typename valueType<T>::type >::type
max( const T & arg ) {
	typedef typename valueType<T>::type V;
	return tools::reduce( arg, [&] (const V& a, const V& b) { return std::max(a, b); }, arg[0] ); 
}

template< class T >
typename activate<T, typename valueType<T>::type >::type
min( const T & arg ) {
	typedef typename valueType<T>::type V;
	return tools::reduce( arg, [&] (const V& a, const V& b) { return std::min(a, b); }, arg[0] ); 
}

/*
 * Comparison operations
 */

template< class T >
inline
typename activate<T, bool>::type
operator==( T const& lhs, T const& rhs )
{
	bool resp = true;
    for (size_t i = 0 ; i < size_impl<T>::size() ; ++i)
        if (lhs[i] != rhs[i]) { resp = false; break; }
    return resp;
}

template< class T >
inline
typename activate<T, bool>::type
operator!=( T const& lhs, T const& rhs )
{
	return (!(lhs == rhs));
}


/*
 * Output operations
 */
template< class T >
typename activate<T, std::ostream &>::type
operator<< (std::ostream &stream, const T & rhs) {
        stream << "[ ";
        for (size_t i = 0; i < size_impl<T>::size() - 1; ++i) 
                stream << rhs[i] << " ; ";
        stream << rhs[size_impl<T>::size()-1] << " ]";
        return stream;
}

// #endif // FLOE_ARITHMETIC_CONTAINER_OPERATORS_HPP

/* vim: set ts=4 sw=4 expandtab : */
