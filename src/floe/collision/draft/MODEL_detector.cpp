// Tentative de modèle pour la détection de collision
// Principe : pouvoir enchevetrer les détecteurs pour simuler les différents niveaux de détection
// La liste des collisions doit être itérable (lazy evaluation)

// GROSSE QST : ça marche comment ce truc ? Si ça ne renvoie que les objets qui sont collisions, comment on détecte les collisions des sous-objets internes à ces objets ?!!
//  * et comment on mixe 2 objets en collisions dans un détecteur sur l'union des sous-objets ?!!
//  * Justement, on ne veut pas gérer les collisions internes aux objets => ça permet de gérer correctement les différents niveaux de détection (par exemple sur différentes unités de calcul)
//  * Par contre, on veut détecter les collisions inter-objets parce que cela demande une communication inter-objets

// QST : comment on spécifie une tolérance ?
//    -> peut-être au dernier niveau de la détection, quand l'algo calcul la distance entre les objets finaux.
//    on peut par exemple fournir une autre fonction distance, prenant en compte une tolérance ?

// RMQ : pour gérer le temps, on pourraît imaginer intercaler un manager d'objet auquel on demanderait de se placer à un instant t donné

// QST : comment on gère la modification des objets ?

// QST : il faut pouvoir récupérer un pas de temps pour prédire les prochaines collisions (notion de vitesse, accélération, rotation, etc ...)
//      -> dans une autre classe non ?
//      un algo naif serait d'avancer d'un pas de temps maximal dt (déterminé en fonction de la vitesse max et de la tolérance) et de
//      faire un genre de dichotomie en cas de détection de collision.
//  => SOL : on peut imaginer 2 étapes de détection de collisions :
//          1) celle avec la meilleure bounding box pour réellement détecter les collisions (on en profite pour récupérer la distance évidemment)
//          2) celle avec une bounding box élargie en fonction de la vitesse de l'objet et le pas de temps maximal afin de détecter les collisions
//              à envisager dans cette échéance. On en profite également pour calculer la distance inter-BB (donc approximative) pour ces "collisions"
//      On aura ainsi une liste couple d'objets et leur distance respective avec une précision adaptée.

// RMQ : la méthode get_filtered_collisions pourrait servir à implémenter le fait qu'on sait que certains objets ne vont pas rentrer en collisions tout de suite
//      => QST : est-ce que le filtre agit par object ou pour chaque couple d'objet ?!!!!

// QST : comment on récupère les points de contacts ?!
//      => dans CollisionList probablement mais il faut imaginer que les fonctions distances puissent renvoyer le point ?

// QST : comment on implémente les détecteurs qui ne sont pas sûrs de leur coup ?
//      du genre, la détection par projection sur droite ...
//      => on pourrait dire que chaque object est entouré par une forme (ici un polygône) équivalente à la capacité de détection de l'algo ?
//          c'est pas trop lourd ?
//      => en tout cas, ça implique d'avoir deux détecteurs au même niveau (plus ou moins)

// RMQ : la détection à l'échelle grossière des collisions ne nécessite pas une précision monstrueuse, on peut donc utiliser des float plutôt que des double (comme dans l'implémentation Sweep&Prune de Danier Joseph Tracy.

template < class Object >
class DetectorModel {
    public :
    DetectorModel();
    ~DetectorModel();

    // Add objects in the detector scope
    template <class Container>
    DetectorModel & append( const Container & objects );

    // Remove objects from the detector scope
    // QST : comment on spécifie les objets à enlever ?
    //  => fonction de filtrage ( => O(n) ) ou sa position dans liste (comment savoir ?) ( => O(1) )
    // QST : do we verify if there is any doubles ?


    /*
    DetectorModel & append( const Object & object );
    DetectorModel & append( Iterator_begin ..., Iterator_end ...);
    */

    // Iterable collision list
    // In another file mayby ? or as children of a base class declared separately ?
    /*
    class CollisionList {
        begin;
        end;
    };
    */
    // At this time, we will use std::vector instead of lazy evaluation
    using CollisionList = std::vector<Object>;

    // Return the list of the detected collisions 
    CollisionList get_collisions();

    // Return the list of the detected collisions in a filtered set of object
    // QST : filter each object of each couple of objects ?
    //         => if we can filter each couple, we can manage collisions between many groups of objects (like parts of two floes)
    // QST : how does the filter function recognizes the object ? unique identifier ? pointer comparison ? => anyway in == definition.
    // QST : filter defined at construction or when collisions detection is asked ?
    CollisionList get_filtered_collisions(const Function & filter);


};

// Une implémentation naive (mais nécessaire à petite échelle) serait de tester la distance entre chaque objets
