/*!
 * \file floe/collision/draft/collision_detector.hpp
 * \brief Draft base class for a detector.
 * \author Roland Denis
 */

#ifndef COLLISION_DETECTOR_HPP_INCLUDED
#define COLLISION_DETECTOR_HPP_INCLUDED

#include <vector>
#include <utility>

namespace floe { namespace collision { namespace draft {

/*! Abstract class for filtering collision.
 * \tparam Object Type of object involved in collision.
 */
template < class Object >
class Filter {
    public :
        /*! Say if the collision must be accepted.
         * \param lhs One object.
         * \param rhs The other object.
         * \return true if the collision is accepted, false otherwise.
         */
        virtual bool pass(const Object & lhs, const Object & rhs) = 0;
};

/*! Filter that accept every collision.
 * \tparam Object Type of object involved in collision.
 */
template < class Object >
class NoFilter : public Filter<Object> {
public :
    //! \return always true ! 
    bool pass(const Object &, const Object &) {
        return true;
    }
};

/*! Base class for a detector of collisions.
 *
 * \tparam Object Type of objects for which collisions will be detected.
 *
 * \todo Possibility to remove objects
 * \remarks No need to specify that objects have move.
 */
template < class Object >
class Detector {
    
public :

    //! Destructor
    virtual ~Detector() = 0;

    /*! Add an object into the detector scope.
     * \param  object The object to add.
     * \return the current instance of the detector.
     */
    virtual Detector<Object> & push_back( const Object & object ) = 0;

    /*! Add a list of objects into the detector scope.
     * \tparam ObjectList   Type of the list, must be iterable.
     * \param  objects        The list of objects to include.
     * \return the current instance of the detector.
     */
    template < class ObjectList >
    Detector<Object> & append( const ObjectList & objects)
    {
        for ( const auto & object : objects ) {
            push_back( object );
        }
        return *this;
    }

    /*! Type of the list of collisions
     * \remark We could later imagine to return a "lazy" list (calculate when the iterator is derefenced)
     */
    using CollisionList = std::vector< std::pair<Object,Object> >;
    
    /*! Return the couples of objects that are in collision.
     * \param filter    A filter that selects which collisions couple are valid.
     */
    virtual CollisionList get_collisions(Filter<Object> & filter = NoFilter<Object>() ) = 0;

    //virtual
    //CollisionList get_collisions(const Filter * filter) = 0;

};

}}} // namespace floe::collision::draft

#endif
