/*!
 * \file floe/state/space_time_state.hpp
 * \brief Space-time state for a floe.
 * \author Roland Denis
 *
 * \todo time state ...
 */

#ifndef FLOE_STATE_SPACE_TIME_STATE_HPP
#define FLOE_STATE_SPACE_TIME_STATE_HPP

#include <cmath>
#include <ostream>

namespace floe { namespace state
{

/*! space-time state of an object
 *
 * It containts no only position and orientation but also speed and rotational
 * speed.
 * Some operations are defined, essentialy to be able to interpolated states
 * and for compatibility with boost::numeric::odeint.
 *
 * \tparam TPos     Position type
 * \tparam TTheta   Orientation type
 * \tparam TSpeed   Speed type
 * \tparam TRot     Rotation type
 *
 */
template <
    typename TPos,
    typename TTheta,
    typename TSpeed = TPos,
    typename TRot = TTheta
>
struct SpaceTimeState
{
    // Type traits
    typedef TPos    Pos;
    typedef TTheta  Theta;
    typedef TSpeed  Speed;
    typedef TRot    Rot;
    
    // State data
    TPos    pos;
    TTheta  theta;
    TSpeed  speed;
    TRot    rot;

    TPos    trans; // total translation (periodic boundaries gap sum)

    //! Compound addition with other state
    template < typename TOtherPos, typename TOtherTheta, typename TOtherSpeed, typename TOtherRot >
    inline
    SpaceTimeState<Pos, Theta, Speed, Rot> & 
    operator += ( SpaceTimeState<TOtherPos, TOtherTheta, TOtherSpeed, TOtherRot> const& other )
    {
        pos     += other.pos;
        theta   += other.theta;
        speed   += other.speed;
        rot     += other.rot;
        return *this;
    }

    //! Compound substraction with other state
    template < typename TOtherPos, typename TOtherTheta, typename TOtherSpeed, typename TOtherRot >
    inline
    SpaceTimeState<Pos, Theta, Speed, Rot> & 
    operator -= ( SpaceTimeState<TOtherPos, TOtherTheta, TOtherSpeed, TOtherRot> const& other )
    {
        pos     -= other.pos;
        theta   -= other.theta;
        speed   -= other.speed;
        rot     -= other.rot;
        return *this;
    }
    
    //! Compound division with other state
    template < typename TOtherPos, typename TOtherTheta, typename TOtherSpeed, typename TOtherRot >
    inline
    SpaceTimeState<Pos, Theta, Speed, Rot> & 
    operator /= ( SpaceTimeState<TOtherPos, TOtherTheta, TOtherSpeed, TOtherRot> const& other )
    {
        pos     /= other.pos;
        theta   /= other.theta;
        speed   /= other.speed;
        rot     /= other.rot;
        return *this;
    }
    
    //! Compound multiplication by a factor
    template < typename T >
    inline
    SpaceTimeState<Pos, Theta, Speed, Rot> &
    operator *= ( T const& factor )
    {
        pos     *= factor;
        theta   *= factor;
        speed   *= factor;
        rot     *= factor;
        return *this;
    }

    //! Compound division by a factor
    template < typename T >
    inline
    SpaceTimeState<Pos, Theta, Speed, Rot> &
    operator /= ( T const& factor )
    {
        pos     /= factor;
        theta   /= factor;
        speed   /= factor;
        rot     /= factor;
        return *this;
    }


    //! Comparison with other state
    template < typename TOtherPos, typename TOtherTheta, typename TOtherSpeed, typename TOtherRot >
    inline
    bool operator == ( SpaceTimeState<TOtherPos, TOtherTheta, TOtherSpeed, TOtherRot> const& other ) const
    {
        return (pos == other.pos) && (theta == other.theta) && (speed == other.speed) && (rot == other.rot);
    }
};


// Deduced operators

//! Addition of states
template < 
    typename LPos, typename LTheta, typename LSpeed, typename LRot, 
    typename RPos, typename RTheta, typename RSpeed, typename RRot 
>
inline
SpaceTimeState<LPos, LTheta, LSpeed, LRot>
operator + ( SpaceTimeState<LPos, LTheta, LSpeed, LRot> lhs, SpaceTimeState<RPos, RTheta, RSpeed, RRot> const& rhs )
{
    lhs += rhs;
    return lhs;
}

//! Substraction of states
template < 
    typename LPos, typename LTheta, typename LSpeed, typename LRot, 
    typename RPos, typename RTheta, typename RSpeed, typename RRot 
>
inline
SpaceTimeState<LPos, LTheta, LSpeed, LRot>
operator - ( SpaceTimeState<LPos, LTheta, LSpeed, LRot> lhs, SpaceTimeState<RPos, RTheta, RSpeed, RRot> const& rhs )
{
    lhs -= rhs;
    return lhs;
}

//! Division of states
template < 
    typename LPos, typename LTheta, typename LSpeed, typename LRot, 
    typename RPos, typename RTheta, typename RSpeed, typename RRot 
>
inline
SpaceTimeState<LPos, LTheta, LSpeed, LRot>
operator / ( SpaceTimeState<LPos, LTheta, LSpeed, LRot> lhs, SpaceTimeState<RPos, RTheta, RSpeed, RRot> const& rhs )
{
    lhs /= rhs;
    return lhs;
}

//! Multiplication by a factor
template < typename TPos, typename TTheta, typename TSpeed, typename TRot, typename T >
inline
SpaceTimeState<TPos, TTheta, TSpeed, TRot>
operator * ( SpaceTimeState<TPos, TTheta, TSpeed, TRot> state, T const& factor )
{
    state *= factor;
    return state;
}

template < typename TPos, typename TTheta, typename TSpeed, typename TRot, typename T >
inline
SpaceTimeState<TPos, TTheta, TSpeed, TRot>
operator * ( T const& factor, SpaceTimeState<TPos, TTheta, TSpeed, TRot> state )
{
    state *= factor;
    return state;
}

//! Division by a factor
template < typename TPos, typename TTheta, typename TSpeed, typename TRot, typename T >
inline
SpaceTimeState<TPos, TTheta, TSpeed, TRot>
operator / ( SpaceTimeState<TPos, TTheta, TSpeed, TRot> state, T const& factor )
{
    state /= factor;
    return state;
}

//! Element-wise absolute value
template < typename TPos, typename TTheta, typename TSpeed, typename TRot>
inline
SpaceTimeState<TPos, TTheta, TSpeed, TRot>
abs ( SpaceTimeState<TPos, TTheta, TSpeed, TRot> state )
{
    state.pos     = std::abs(state.pos);
    state.theta   = std::abs(state.theta);
    state.speed   = std::abs(state.speed);
    state.rot     = std::abs(state.rot);
    return state;
}

//! Stream state
template < typename TPos, typename TTheta, typename TSpeed, typename TRot>
inline
std::ostream& operator<< (std::ostream & out, SpaceTimeState<TPos, TTheta, TSpeed, TRot> const& state)
{
    out          << "Position: " << state.pos
        << " ; " << "Orientation: " << state.theta
        << " ; " << "Speed: " << state.speed
        << " ; " << "Rotation: " << state.rot;
    return out;
}


}} // namespace floe::state

namespace std
{
    using floe::state::abs;
} // namespace std
#endif // FLOE_STATE_SPACE_TIME_STATE_HPP

