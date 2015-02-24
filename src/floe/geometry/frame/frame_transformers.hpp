/*!
 * \file floe/geometry/frame/frame_transformers.hpp
 * \brief Transformation (translate, rotate, scale, ...) for frames.
 * \author Roland Denis
 * \see <a href="http://www.boost.org/doc/libs/1_57_0/libs/geometry/doc/html/geometry/reference/algorithms/transform/transform_3_with_strategy.html">boost/geometry/algorithms/transform.hpp</a>
 */

#ifndef FLOE_GEOMETRY_FRAME_FRAME_TRANSFORMERS_HPP
#define FLOE_GEOMETRY_FRAME_FRAME_TRANSFORMERS_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/geometry/strategies/transform/matrix_transformers.hpp>
#include <boost/geometry/util/select_most_precise.hpp>

#include "floe/geometry/core/access.hpp"

namespace floe { namespace geometry { namespace frame
{


// Some tools
namespace {

namespace bg = boost::geometry::strategy::transform;

//! Matrix type
template < typename T >
using matrix_type = boost::numeric::ublas::matrix<T>;

//! Translate operation
template <typename T>
inline
matrix_type<T> translate( T const& x, T const& y )
{
    return bg::translate_transformer<T,2,2>(x,y).matrix();
}

/*! Rotation operator
 *
 * The boost rotate_transformer is clockwise oriented.
 */
template <typename T>
inline
matrix_type<T> rotate( T const& theta )
{
    return bg::rotate_transformer<boost::geometry::radian,T,2,2>(-theta).matrix();
}

} // namespace

//! Transformer type
template < typename T >
using transformer_type = bg::ublas_transformer<T,2,2>;

/*! Transformer strategy from a frame to the canonical frame
 *
 * \tparam  TFrame  Type of the frame (auto-deducted)
 * \param   frame   The frame.
 * \return  transfomer whose matrix equals translate(center)*rotate(theta)
 */
template < 
    typename TFrame,
    typename T = typename TFrame::coordinate_type
>
inline
transformer_type<T> transformer( TFrame const& frame )
{
    using floe::geometry::get;
    return transformer_type<T>(
        prod(
            translate(get<0>(frame.center()), get<1>(frame.center())),
            rotate(frame.theta())
        )
    );
}

/*! Transformer strategy from the canonical frame to a given frame
 *
 * \tparam  TFrame  Type of the frame (auto-deducted).
 * \param   frame   The frame.
 * \return  transformer whose matrix equals rotate(-theta)*translate(-center)
 */
template < 
    typename TFrame,
    typename T = typename TFrame::coordinate_type
>
inline
transformer_type<T> itransformer( TFrame const& frame )
{
    using floe::geometry::get;
    return transformer_type<T>(
        prod( 
            rotate(-frame.theta()),
            translate(-get<0>(frame.center()), -get<1>(frame.center()))
        )
    );
}

/*! Transformer strategy for a frame to an another
 *
 * \tparam  TFrame1 Type of the source frame.
 * \tparam  TFrame2 Type of the destination frame.
 * \param   frame1  The source frame.
 * \param   frame2  The destination frame.
 *
 * \return transformer whose matrix equals rotate(-theta2)*translate(center1-center2)*rotate(theta1)
 */
template < 
    typename TFrame1, 
    typename TFrame2,
    typename T1 = typename TFrame1::coordinate_type,
    typename T2 = typename TFrame2::coordinate_type,
    typename TOutput = typename boost::geometry::select_most_precise<T1,T2>::type
>
inline
transformer_type<TOutput> transformer( TFrame1 const& frame1, TFrame2 const& frame2 )
{
    using floe::geometry::get;
    return transformer_type<TOutput>(
        prod(
            rotate(-frame2.theta()),
            matrix_type<TOutput>(prod(
                translate( 
                    get<0>(frame1.center()) - get<0>(frame2.center()),
                    get<1>(frame1.center()) - get<1>(frame2.center())
                ),
                rotate(frame1.theta())
            ))
        )
    );
}

/*! Transformer strategy for a frame to an another frame for objects in canonical frame
 *
 * \tparam  TFrame1 Type of the source frame.
 * \tparam  TFrame2 Type of the destination frame.
 * \param   frame1  The source frame.
 * \param   frame2  The destination frame.
 *
 * \return transformer whose matrix equals translate(center2)*rotate(theta2-theta1)*translate(-center1)
 */
template < 
    typename TFrame1, 
    typename TFrame2,
    typename T1 = typename TFrame1::coordinate_type,
    typename T2 = typename TFrame2::coordinate_type,
    typename TOutput = typename boost::geometry::select_most_precise<T1,T2>::type
>
inline
transformer_type<TOutput> itransformer( TFrame1 const& frame1, TFrame2 const& frame2 )
{
    using floe::geometry::get;
    return transformer_type<TOutput>(
        prod(
            translate( get<0>(frame2.center()), get<1>(frame2.center()) ),
            matrix_type<TOutput>(prod(
                rotate(frame2.theta()-frame1.theta()),
                translate( -get<0>(frame1.center()), -get<1>(frame1.center()) )
            ))
        )
    );
}
}}} // namespace floe::geometry::frame

#endif // FLOE_GEOMETRY_FRAME_FRAME_TRANSFORMERS_HPP

