/*!
 * \file floe/io/matlab/list_so.hpp
 * \brief Structures reproducing the list_so from the Matlab code. Only for I/O purpose.
 * \author Roland Denis
 *
 * The name of the structures are (more or less) the same as in the matlab code.
 */

#ifndef FLOE_IO_MATLAB_LIST_SO_HPP
#define FLOE_IO_MATLAB_LIST_SO_HPP

#include <vector>
#include <cstddef>
#include <array>
#include <boost/multi_array.hpp>

namespace floe { namespace io { namespace matlab
{

// Typedefs and structures

template <typename T> using point_type = std::array<T, 2>;
template <typename T> using mark_type = std::array<point_type<T>, 3>;
template <typename T> using ssdif_type = point_type<T>;

template <typename T>
struct border_type {
    std::vector< point_type<T> > S;
    std::vector< ssdif_type<T> > ssdif;
};

template <typename T>
struct mesh_type {
    std::vector< point_type<T> >     vertices;
    std::vector< point_type<T> >     border;
    std::vector< std::array<T,3> > triangles;
};

template <typename T>
struct own_disk_type {
    point_type<T> center_disk;
};

template <typename T>
struct pt_own_mark_type {
    mesh_type<T>       S;
    own_disk_type<T>   D;
};

template <typename T>
struct disk_rep_type {
    T  id;
    T  radius;
};

template <typename T>
struct disks_type {
    T nC;
    std::vector< disk_rep_type<T> > rep;
};

template <typename T>
struct quadrature_type {
    std::array<point_type<T>,3> points;
};

/*! Structure reproducing Solid structure
 *
 * \tparam T Value type used in Matlab
 */
template <
    typename T = double
>
struct MatlabSolid
{
    // Typedefs and structures
    typedef T value_type;

    // Datas
    value_type              tau;
    value_type              size_maille;
    point_type<T>           center_mass;
    mark_type<T>            mark;
    border_type<T>          border;
    pt_own_mark_type<T>     pt_own_mark;
    point_type<T>           center_disk;
    value_type              radius;
    value_type              area;
    disks_type<T>           disks;
    value_type          d_c;
};

/*! Structure reproducing Movement structure
 *
 * \tparam T Value type used in Matlab
 */
template <
    typename T = double
>
struct MatlabMovement
{
    // Typedefs and structures
    typedef T value_type;

    // Datas
    point_type<T>      center;
    point_type<T>      vcenter;
    value_type      theta;
    value_type      vtheta;
    std::vector< quadrature_type<T> >    pt_reel;
    std::vector<T>     detJ;
    value_type      mass_tot;
    std::vector<T>     mass;
    std::vector<T>     area;
    value_type      denom;
    value_type      time;
    value_type      timecurr;
    value_type      obs;
};


/*! Class reproducing List_Solid structure
 *
 * \tparam T Value type used in Matlab
 */
template <
    typename T = double
>
struct MatlabListSolid
{
    // Typedefs and structures
    typedef T value_type;
    typedef boost::multi_array<value_type,2> array_type;
    typedef std::array<value_type,3> contact_type; // ????
    typedef std::vector<contact_type>  contact_list_type;
    typedef boost::multi_array<contact_list_type,2> contact_array_type;


    // Datas
    array_type  delta_t;
    array_type  indic;
    array_type  d_secu;
    array_type  d_opt;
    contact_array_type point_contact;
    std::vector< MatlabSolid<value_type> > geo;
    std::vector< MatlabMovement<value_type> > mov;
    value_type cinetic;

    // Default constructor
    MatlabListSolid()
        :   delta_t{ boost::extents[0][0] },
            indic{ boost::extents[0][0] },
            d_secu{ boost::extents[0][0] },
            d_opt{ boost::extents[0][0] },
            point_contact{ boost::extents[0][0] }
    {}
};


}}} // namespace floe::io::matlab

#endif // FLOE_MATLAB_LIST_SO_HPP
