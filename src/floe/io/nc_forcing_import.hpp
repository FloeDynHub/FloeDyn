/*!
 * \file floe/io/nc_forcing_import.hpp
 * \brief Import of inhomogeneous ocean/wind forcing from a NetCDF4 file (HDF5 backend).
 * \author Quentin Jouet
 */

#ifndef FLOE_IO_NC_FORCING_IMPORT_HPP
#define FLOE_IO_NC_FORCING_IMPORT_HPP

#include <string>
#include <vector>
#include <cstdint>
#include "H5Cpp.h"

namespace floe { namespace io {

struct NcForcingData {
    std::vector<double>  time;    // seconds, shape (nt,)
    std::vector<int64_t> x;       // grid indices, shape (nx,)
    std::vector<int64_t> y;       // grid indices, shape (ny,)
    double delta_x;               // meters per x grid unit
    double delta_y;               // meters per y grid unit
    std::size_t nt, ny, nx;
    // flattened [t * ny * nx + iy * nx + ix]
    std::vector<double> wind_x, wind_y;
    std::vector<double> ocean_x, ocean_y;
};

inline NcForcingData read_nc_forcing(std::string const& filename)
{
    using namespace H5;
    NcForcingData d;
    H5File file(filename, H5F_ACC_RDONLY);

    auto read_1d_double = [&](std::string const& name, std::vector<double>& out, std::size_t& n) {
        DataSet ds = file.openDataSet(name);
        hsize_t dims[1];
        ds.getSpace().getSimpleExtentDims(dims);
        n = dims[0];
        out.resize(n);
        ds.read(out.data(), PredType::NATIVE_DOUBLE);
        return ds;
    };

    auto read_1d_int64 = [&](std::string const& name, std::vector<int64_t>& out, std::size_t& n) {
        DataSet ds = file.openDataSet(name);
        hsize_t dims[1];
        ds.getSpace().getSimpleExtentDims(dims);
        n = dims[0];
        out.resize(n);
        ds.read(out.data(), PredType::NATIVE_INT64);
        return ds;
    };

    auto read_str_attr = [](DataSet const& ds, std::string const& attr_name) {
        Attribute attr = ds.openAttribute(attr_name);
        std::string val;
        attr.read(attr.getStrType(), val);
        return val;
    };

    read_1d_double("time", d.time, d.nt);

    {
        auto ds = read_1d_int64("x", d.x, d.nx);
        d.delta_x = std::stod(read_str_attr(ds, "delta_x"));
    }
    {
        auto ds = read_1d_int64("y", d.y, d.ny);
        d.delta_y = std::stod(read_str_attr(ds, "delta_y"));
    }

    auto read_3d = [&](std::string const& name) {
        DataSet ds = file.openDataSet(name);
        std::vector<double> out(d.nt * d.ny * d.nx);
        ds.read(out.data(), PredType::NATIVE_DOUBLE);
        return out;
    };

    d.wind_x  = read_3d("wind_xvelocity");
    d.wind_y  = read_3d("wind_yvelocity");
    d.ocean_x = read_3d("ocean_xvelocity");
    d.ocean_y = read_3d("ocean_yvelocity");

    std::cout << "Loaded NetCDF forcing: " << d.nt << " time steps, grid "
              << d.ny << "x" << d.nx
              << " (delta_x=" << d.delta_x << "m, delta_y=" << d.delta_y << "m)" << std::endl;
    return d;
}

}} // namespace floe::io

#endif // FLOE_IO_NC_FORCING_IMPORT_HPP
