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


    template <class G>
    struct py_grid_signature
    {
        // using neighbors_count = typename G::neighbors_count_type& (G::*)(const typename
        // G::size_type&) const;
        using neighbors_indices
            = typename G::neighbors_indices_type (G::*)(const typename G::size_type&);
        // using neighbors_distances = typename G::neighbors_distances_type& (G::*)(const typename
        // G::size_type&) const;
    };
}


#endif  // PYFASTSCAPELIB_GRID_H
