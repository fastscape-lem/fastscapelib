#ifndef PYFASTSCAPELIB_GRID_H
#define PYFASTSCAPELIB_GRID_H

#include <functional>
#include <stdexcept>

#include "pybind11/pybind11.h"

#include "fastscapelib/grid.hpp"
#include "fastscapelib/profile_grid.hpp"
#include "fastscapelib/raster_grid.hpp"
#include "fastscapelib/unstructured_mesh.hpp"

#include "pytensor_utils.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


namespace fastscapelib
{

    using py_profile_grid = fs::profile_grid_xt<fs::pyarray_selector>;
    using py_raster_grid = fs::raster_grid_xt<fs::pyarray_selector, fs::raster_connect::queen>;
    using py_unstructured_mesh = fs::unstructured_mesh_xt<fs::pyarray_selector>;

    template <class G>
    struct py_grid_funcs
    {
    public:
        using size_type = typename G::size_type;
        using count_type = typename G::neighbors_count_type;
        using indices_type = typename G::neighbors_indices_type;
        using distances_type = typename G::neighbors_distances_type;
        using neighbors_type = typename G::neighbors_type;

        py_grid_funcs() = default;

        std::function<count_type(const G&, size_type)> m_neighbors_count
            = [this](const G& g, size_type idx)
        {
            check_in_bounds(g, idx);
            return g.neighbors_count(idx);
        };

        // no const since it may update the grid cache internally
        std::function<indices_type(G&, size_type)> m_neighbors_indices = [this](G& g, size_type idx)
        {
            check_in_bounds(g, idx);
            return g.neighbors_indices(idx);
        };

        std::function<distances_type(const G&, size_type)> m_neighbors_distances
            = [this](const G& g, size_type idx)
        {
            check_in_bounds(g, idx);
            return g.neighbors_distances(idx);
        };

        // no const since it may update the grid cache internally
        std::function<neighbors_type(G&, size_type)> m_neighbors = [this](G& g, size_type idx)
        {
            check_in_bounds(g, idx);
            return g.neighbors(idx);
        };

    private:
        inline void check_in_bounds(const G& g, size_type idx)
        {
            if (idx < 0 || idx >= g.size())
            {
                throw std::out_of_range("grid index out of range");
            }
        }
    };

    template <class G>
    void add_neighbor_methods(py::class_<G>& pyg)
    {
        auto grid_funcs = py_grid_funcs<G>();

        pyg.def("neighbors_count", grid_funcs.m_neighbors_count)
            .def("neighbors_indices", grid_funcs.m_neighbors_indices)
            .def("neighbors_distances", grid_funcs.m_neighbors_distances)
            .def("neighbors", grid_funcs.m_neighbors);
    }

}

#endif  // PYFASTSCAPELIB_GRID_H
