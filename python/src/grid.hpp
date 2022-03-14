#ifndef PYFASTSCAPELIB_GRID_H
#define PYFASTSCAPELIB_GRID_H

#include <functional>

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
    struct py_grid_funcs
    {
        using size_type = typename G::size_type;
        using count_type = typename G::neighbors_count_type;
        using indices_type = typename G::neighbors_indices_type;
        using distances_type = typename G::neighbors_distances_type;
        using neighbors_type = typename G::neighbors_type;

        std::function<count_type(const G&, size_type)> neighbors_count
            = [](const G& g, size_type idx) { return g.neighbors_count(idx); };

        // no const since it may update the grid cache internally
        std::function<indices_type(G&, size_type)> neighbors_indices
            = [](G& g, size_type idx) { return g.neighbors_indices(idx); };

        std::function<distances_type(const G&, size_type)> neighbors_distances
            = [](const G& g, size_type idx) { return g.neighbors_distances(idx); };

        // no const since it may update the grid cache internally
        std::function<neighbors_type(G&, size_type)> neighbors
            = [](G& g, size_type idx) { return g.neighbors(idx); };
    };
}

#endif  // PYFASTSCAPELIB_GRID_H
