// Tentative de modèle pour gérer les niveaux de détecteurs

#include "collision/detector.hpp"

template<
    ObjectList,
    Object
>
detector_manager( const ObjectList & objects ) {

    Detector<Object> dectector;
    detector.append(objects);

    ...

    detector.update();
    auto col_list = detector.get_collisions();
    std::vector<Object> tmp_objects;

    for ( auto & col : col_list ) {
        auto sub_objects = get_intersection( col.first, col.second );
        if ( sub_objects.size() > 0 ) {
            tmp_objects.push_back({sub_objects});
        }
    }

    // Manage next level collisions in each object of objects and tmp_objects
    ...
}
