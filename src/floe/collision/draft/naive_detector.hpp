/*!
 * \file floe/collision/draft/naive_detector.hpp
 * \brief Naive implementation of a detector based on the draft base class.
 * \author Roland Denis
 */

#ifndef COLLISION_NAIVE_DETECTOR_HPP_INCLUDED
#define COLLISION_NAIVE_DETECTOR_HPP_INCLUDED

#include <vector>

#include "collision/detector.h"

namespace floe { collision { draft {

/*! Naive (but usefull) implementation of a detector.
 *
 * It tests collision for every couple of objects ( O(n²) ).
 *
 * \tparam Object Type of object involved in collision.
 *
 * \see collision::Detector
 */
template < class Object >
class NaiveDetector : public Detector<Object> 
{

private :
    std::vector<const Object> m_objects; //! Internal list of objects

public :
    
    //! Constructor
    NaiveDetector() {}
    
    //! Destructor
    ~NaiveDetector() {}
    
    /*! Add an object into the detector scope.
     * \param object The object to add.
     * \return the current instance of the detector.
     */
    Detector<Object> & push_back( const Object & object )
    {
        m_objects.push_back(object);
        return *this;
    }

    //! Type of the list of collisons.
    using CollisionList = typename Detector<Object>::CollisionList;

    /*! Return the couples of objects that are in collision.
     * \param filter    A filter that selects which collisions couple are valid.
     */
    CollisionList get_collisions(Filter<Object> & filter = NoFilter<Object>() ) {
        CollisionList col_list;

        using std::begin;
        using std::end;
        using std::next;
        using std::prev;

        // For every couple ...
        // Remark: we could use constant iterator instead ...
        for ( auto lhs_it = begin(m_objects); lhs_it != prev(end(m_objects)); ++lhs_it ) {
            for ( auto rhs_it = next(lhs_it); rhs_it != end(m_objects); ++rhs_it ) {

                if ( filter->pass( *lhs_it, *rhs_it) && collide( *lhs_it, *rhs_it ) ) {
                    col_list.push_back( { *lhs_it, *rhs_it } );
                }
            }
        }

        return col_list;
    }


};

}}} // namespace floe::collision::draft

#endif
