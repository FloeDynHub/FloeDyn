#ifndef FLOE_GEOMETRY_GEOMETRIES_MULTI_SSP_POINT_CLOUD
#define FLOE_GEOMETRY_GEOMETRIES_MULTI_SSP_POINT_CLOUD

#include <cstddef>
#include <type_traits>

#include <boost/concept/assert.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "floe/geometry/geometries/concepts/multi_point_concept.hpp"
#include "floe/geometry/geometries/concepts/multi_simple_static_polygon_concept.hpp"
#include "floe/geometry/core/static_num_points.hpp"

namespace floe { namespace geometry 
{

namespace {
    template <typename T>
    inline T* as_ptr(T* obj) { return obj; }
    
    template <typename T>
    inline T const* as_ptr(T const* obj) { return obj; }

    template <typename T>
    inline T* as_ptr(T& obj) { return &obj; }

    template <typename T>
    inline T const* as_ptr(T const& obj) { return &obj; }

} // namespace

/*! Constant Multi Simple Static Polygon from a point cloud
 *
 * Given a multi-point and a connectivity list, it models
 * a ConstMultiSimpleStaticPolygon concept. 
 * The ConstSimpleStaticPolygon model that will be returned through
 * the range concept must be specified by template.
 *
 * Point list and connectivity list can be passed by const-ref by
 * specifying the right template.
 *
 * \tparam TPolygon         Type of the ConstSimpleStaticPolygon model.
 * \tparam TMultiPoint      Type of a MultiPoint model. 
 * \tparam TConnectivity    Type of the connectivity list.
 */
template <
    typename TPolygon,
    typename TMultiPoint,
    typename TConnectivity
>
class MultiSSPPointCloud
{

    // Concept checking
    BOOST_CONCEPT_ASSERT( (concepts::ConstSimpleStaticPolygon< typename std::remove_pointer<TPolygon>::type >) );
    BOOST_CONCEPT_ASSERT( (concepts::ConstMultiPoint< typename std::remove_pointer<TMultiPoint>::type >) );
    BOOST_CONCEPT_ASSERT( (boost::RandomAccessRangeConcept< typename std::remove_pointer<TConnectivity>::type >) );
    // TODO: Concept for TConnectivity::value_type

    // type traits
    typedef typename std::remove_pointer<typename std::decay<TMultiPoint>::type>::type multi_point_type;
    typedef typename std::remove_pointer<typename std::decay<TConnectivity>::type>::type connectivity_type;

public:
    template< typename TMP, typename TC >
    MultiSSPPointCloud( TMP && points, TC && connect )
        : m_points{std::forward<TMP>(points)}, m_connect{std::forward<TC>(connect)}
    {}

    inline std::size_t size() const  { return m_connect->size(); }

private:
    class PolygonBuilder 
    {
        public:
            PolygonBuilder() {}
            PolygonBuilder(multi_point_type const* cloud) : m_cloud(cloud) {}
            TPolygon operator() (typename connectivity_type::value_type const& indexes) const;
        private:
            multi_point_type const* m_cloud;
    };

public:

    ////// Type Traits //////

    //! Const iterator type
    typedef boost::transform_iterator<
            PolygonBuilder,
            typename connectivity_type::const_iterator>
        const_iterator;

    //! Const reverse iterator type
    typedef boost::transform_iterator<
            PolygonBuilder,
            typename connectivity_type::const_reverse_iterator>
        const_reverse_iterator;

    typedef typename const_iterator::value_type         value_type;
    typedef typename connectivity_type::size_type       size_type;
    typedef typename connectivity_type::difference_type difference_type;
    typedef typename const_iterator::reference          reference;
    
    ////// Iterators //////
    //! forward begin
    inline
    const_iterator begin() const 
    {
        return const_iterator{as_ptr(m_connect)->begin(), PolygonBuilder{as_ptr(m_points)}};
    }

    //! forward end
    inline
    const_iterator end() const
    {
        return const_iterator{as_ptr(m_connect)->end(), PolygonBuilder{as_ptr(m_points)}};
    }
    
    //! reverse begin
    inline
    const_iterator rbegin() const
    {
        return const_reverse_iterator{as_ptr(m_connect)->rbegin(), PolygonBuilder{as_ptr(m_points)}};
    }

    //! reverse end
    inline
    const_iterator rend() const
    {
        return const_reverse_iterator{as_ptr(m_connect)->rend(), PolygonBuilder{as_ptr(m_points)}};
    }

    //! random access
    inline
    TPolygon operator[] (std::size_t i) const
    {
        return PolygonBuilder(as_ptr(m_points))((*as_ptr(m_connect))[i]);
    }
    
private:
    TMultiPoint m_points;
    TConnectivity m_connect;

};


/*
 * Some internal tools
 */
namespace {
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

template <
    typename TPolygon,
    typename TMultiPoint,
    typename TIndex,
    std::size_t... I
>
TPolygon build_polygon_impl( TMultiPoint const& cloud, TIndex const& indexes, seq<I...> )
{
    return TPolygon{ cloud[indexes[I]]... };
}

template < 
    typename TPolygon,
    typename TMultiPoint,
    typename TIndex
>
TPolygon build_polygon( TMultiPoint const& cloud, TIndex const& indexes )
{
    return build_polygon_impl<TPolygon>( 
        cloud, indexes, 
        typename gen_seq< static_num_points<TPolygon>::value >::type()
    );
}

} // namespace

template <
    typename TPolygon,
    typename TMultiPoint,
    typename TConnectivity
>
TPolygon
MultiSSPPointCloud<TPolygon,TMultiPoint,TConnectivity>::PolygonBuilder::
operator() ( typename connectivity_type::value_type const& indexes ) const
{
    return build_polygon<TPolygon>( *m_cloud, indexes );
}

}} // namespace floe::geometry

namespace boost { namespace geometry { namespace traits
{

//! Tag
template <
    typename TPolygon,
    typename TMultiPoint,
    typename TConnectivity
>
struct tag< floe::geometry::MultiSSPPointCloud<TPolygon,TMultiPoint,TConnectivity> >
{
    typedef multi_simple_static_polygon_tag type;
};

}}} // namespace boost::geometry::traits

namespace boost
{

template <
    typename TPolygon,
    typename TMultiPoint,
    typename TConnectivity
>
struct range_value< floe::geometry::MultiSSPPointCloud<TPolygon,TMultiPoint,TConnectivity> >
{
    typedef TPolygon type;
};

} // namespace boost


#endif // FLOE_GEOMETRY_GEOMETRIES_MULTI_SSP_POINT_CLOUD
