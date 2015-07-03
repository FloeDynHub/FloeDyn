/*!
 * \file variable/floe_group.hpp
 * \brief Floe Configuration class
 * \author Quentin Jouet
 */

#ifndef VARIABLE_FLOES_HPP
#define VARIABLE_FLOES_HPP

#include <iostream>
#include "floe/variable/floe_group_h.hpp"

#include "floe/io/matlab/list_so_to_floes.hpp"
#include "floe/io/matlab/list_so_import.hpp"
#include "floe/io/matlab/list_so.hpp"

#include "H5Cpp.h"
#include <algorithm>
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif


namespace floe { namespace variable
{

/*! FloeGroup
 *
 * It represents the set of Floes.
 *
 */


template <
    typename TFloe
>
class FloeGroup
{

public:

    using floe_type = TFloe;
    using floe_group_h_type = floe::variable::FloeGroup_h<
        typename floe_type::floe_h_type
    >;
    using value_type = typename TFloe::value_type;

    //! Default constructor.
    FloeGroup() : m_out_file{nullptr} {}


    void load_matlab_config(std::string filename);

    void out_hdf5(value_type time);

    inline void add_floe( floe_type& floe )
        {
            m_list_floe.push_back(std::move(floe));
            // auto& f = m_list_floe.back(); // TODO Why is this different from *1 ?
            // f.update();
            // m_floe_group_h.add_floe(f.get_floe_h());
        }

    inline floe_group_h_type const& get_floe_group_h() const { return m_floe_group_h; }
    inline floe_group_h_type& get_floe_group_h() { return m_floe_group_h; }
    inline std::vector<floe_type>& get_floes() { return m_list_floe; }

    //! kinetic energy of the group
    value_type kinetic_energy();

private:

    std::vector<floe_type> m_list_floe;
    floe_group_h_type m_floe_group_h;
    H5File* m_out_file;

};


template <
    typename TFloe
>
void FloeGroup<TFloe>::load_matlab_config(std::string filename) {
    using namespace floe::io::matlab;
    MatlabListSolid<double> list_so;
    cout << "Reading \"" << filename << "\" ... " << endl;
    read_list_so_from_file( filename, list_so);

    cout << "Importing floes ... " << endl;
    // m_list_floe = list_so_to_floes<floe_type>( list_so );
    for ( auto& floe : list_so_to_floes<floe_type>( list_so ) )
        add_floe(floe);
    for ( auto& floe : m_list_floe ){ // *1 ; todo : avoid that step
        floe.update();
        m_floe_group_h.add_floe(floe.get_floe_h());
    }
};


template <
    typename TFloe
>
void FloeGroup<TFloe>::out_hdf5(value_type time) {
    try
    {   
        const H5std_string  FILE_NAME( "out/out.h5" );
        const int   NY = 2;
        const int   RANK = 2;
        /*
         * Turn off the auto-printing when failure occurs so that we can
         * handle the errors appropriately
         */
        Exception::dontPrint();
        /*
         * Create a new file using H5F_ACC_TRUNC access,
         * default file creation properties, and default file
         * access properties.
         */
        if (m_out_file == nullptr)
            m_out_file = new H5File( FILE_NAME, H5F_ACC_TRUNC );
        H5File& file( *m_out_file );
        /*
         * Define the size of the array and create the data space for fixed
         * size dataset.
         */

        FloatType datatype( PredType::NATIVE_DOUBLE );
        datatype.setOrder( H5T_ORDER_LE );
        hsize_t     dimsf[2];              // dataset dimensions
        dimsf[1] = NY;

        Group group(file.createGroup(H5std_string{std::to_string(time)}));

        for (std::size_t i=0; i!=m_list_floe.size(); ++i)
        {
            auto& boundary = m_list_floe[i].geometry().outer();
            dimsf[0] = boundary.size();
            DataSpace dataspace( RANK, dimsf );
            /*
             * Create a nFew dataset within the file using defined dataspace and
             * datatype and default dataset creation properties.
             */
            DataSet dataset = group.createDataSet(H5std_string{std::to_string(i)},datatype, dataspace);
            value_type data[dimsf[0]][dimsf[1]];
            for (std::size_t j = 0; j!= dimsf[0]; ++j)
            {
                data[j][0] = boundary[j].x;
                data[j][1] = boundary[j].y;
            }
            /*
             * Write the data to the dataset using default memory space, file
             * space, and transfer properties.
             */
            dataset.write( data, PredType::NATIVE_DOUBLE );
        }

    }  // end of try block
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
        error.printError();
        // return -1;
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
        error.printError();
        // return -1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
        error.printError();
        // return -1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException error )
    {
        error.printError();
        // return -1;
    }
};


template<typename TFloe>
typename FloeGroup<TFloe>::value_type
FloeGroup<TFloe>::kinetic_energy()
{
    return std::accumulate(
        m_list_floe.begin(), m_list_floe.end(), 0. , 
        [](value_type partial_sum, floe_type& floe) {return partial_sum + floe.kinetic_energy(); }
    );
}


}} // namespace floe::variable


#endif // VARIABLE_FLOES_HPP
