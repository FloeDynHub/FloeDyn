/*!
 * \file ope/explicit_physical_data.hpp
 * \brief Physical datas (air & ocean speed)
 * \author Quentin Jouet
 */

#ifndef OPE_PHYSICAL_DATA_HPP
#define OPE_PHYSICAL_DATA_HPP

#include "floe/geometry/arithmetic/point_operators.hpp"
#include "floe/io/matlab/topaz_import.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream> // DEBUG


namespace floe { namespace ope
{

/*! PhysicalData
 *
 * Gives air and ocean datas throughout time and space
 *
 */


template <typename TPoint>
class PhysicalData
{

public:

    using point_type = TPoint;
    using value_type = decltype(TPoint::x);
    using point_vector = std::vector<point_type>;

    //! Constructor
    PhysicalData(value_type const& time_ref) : m_time_ref{time_ref} {}

    //! water speed accessor
    // point_type water_speed(point_type pt = {0,0}) { return {0,0}; } // Null
    // point_type water_speed(point_type pt = {0,0}) { return - pt / std::max(2000., norm2(pt)); } // convergent current
    point_type water_speed(point_type pt = {0,0}) {
        int sign{(int(m_time_ref / 3000) % 3 == 2) ? 1 : -1};
        // return sign * pt / (3 * norm2(pt));
        return sign * pt / std::max(600., 3 * norm2(pt));
    } // convergent current

    //! air speed accessor
    point_type air_speed(point_type = {0,0}) { return {0,0}; } // Null
    // point_type air_speed(point_type = {0,0}) { return {10,0}; } // Null

    void update_water_speed(point_type diff_speed){}
    //! OBL speed accessor for output
    inline point_type OBL_speed() const { return {0,0}; }

    //! Load ocean and wind data from a topaz file
    void load_matlab_topaz_data(std::string const& filename) {}

private:

    value_type const& m_time_ref; //!< reference to time variable
};


}} // namespace floe::ope


#endif // OPE_PHYSICAL_DATA_HPP