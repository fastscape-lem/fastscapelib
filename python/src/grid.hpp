#ifndef PYFASTSCAPELIB_GRID_H
#define PYFASTSCAPELIB_GRID_H

#include "fastscapelib/grid.hpp"
#include "fastscapelib/profile_grid.hpp"
#include "fastscapelib/raster_grid.hpp"
#include "fastscapelib/unstructured_mesh.hpp"

#include "pytensor_utils.hpp"


namespace fs = fastscapelib;


namespace fastscapelib
{

    using py_profile_grid = fs::profile_grid_xt<fs::pyarray_selector>;
    using py_raster_grid = fs::raster_grid_xt<fs::pyarray_selector, fs::raster_connect::queen>;
    using py_unstructured_mesh = fs::unstructured_mesh_xt<fs::pyarray_selector>;
}


#endif  // PYFASTSCAPELIB_GRID_H
