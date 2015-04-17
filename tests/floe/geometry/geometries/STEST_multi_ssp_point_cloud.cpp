#include <iostream>
#include <array>
#include <vector>
#include <cstddef>

#include <boost/concept/assert.hpp>

#include "floe/geometry/geometry.hpp"
#include "floe/geometry/geometries/point.hpp"
#include "floe/geometry/geometries/triangle.hpp"
#include "floe/geometry/geometries/multi_ssp_point_cloud.hpp"

#include "floe/geometry/geometries/concepts/multi_simple_static_polygon_concept.hpp"

#include <boost/geometry/multi/geometries/multi_point.hpp>


int main()
{
    using namespace std;

    namespace fg = floe::geometry;
    using real = double;
    using Point = fg::Point<real>;
    using Triangle = fg::Triangle<Point>;
    using MP = boost::geometry::model::multi_point<Point>;
    using Connect = std::vector< std::array<std::size_t,3> >;
    using MSSP_ref = const fg::MultiSSPPointCloud<Triangle, MP const&, Connect const&>;
    using MSSP_copy = const fg::MultiSSPPointCloud<Triangle, MP, Connect>;
    using MSSP_ptr = const fg::MultiSSPPointCloud<Triangle, MP const*, Connect const*>;

    BOOST_CONCEPT_ASSERT( (boost::RandomAccessRangeConcept<MSSP_ref>) );
    BOOST_CONCEPT_ASSERT( (boost::RandomAccessRangeConcept<MSSP_copy>) );
    BOOST_CONCEPT_ASSERT( (boost::RandomAccessRangeConcept<MSSP_ptr>) );

    BOOST_CONCEPT_ASSERT( (fg::concept::ConstMultiSimpleStaticPolygon<MSSP_ref>) );
    BOOST_CONCEPT_ASSERT( (fg::concept::ConstMultiSimpleStaticPolygon<MSSP_copy>) );
    BOOST_CONCEPT_ASSERT( (fg::concept::ConstMultiSimpleStaticPolygon<MSSP_ptr>) );

    MP cloud;
    cloud.push_back({0,0});
    cloud.push_back({0,1});
    cloud.push_back({1,0});
    cloud.push_back({1,1});
    
    Connect connect;
    connect.push_back({{ 0, 1, 2}});
    connect.push_back({{ 1, 2, 3}});
    connect.push_back({{0,2,3}});

    MSSP_ref mesh { cloud, connect };

    for ( const auto& triangle : mesh )
        cout << fg::dsv(triangle) << endl;

    MSSP_copy mesh2 { cloud, connect };

    for ( const auto& triangle : mesh2 )
        cout << fg::dsv(triangle) << endl;

    MSSP_ptr mesh3 { &cloud, &connect };

    for ( const auto& triangle : mesh3 )
        cout << fg::dsv(triangle) << endl;
    return 0;
}
