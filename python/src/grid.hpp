#ifndef PYFASTSCAPELIB_GRID_H
#define PYFASTSCAPELIB_GRID_H

#include <functional>
#include <stdexcept>

#include "pybind11/pybind11.h"

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/grid/profile_grid.hpp"
#include "fastscapelib/grid/raster_grid.hpp"
#include "fastscapelib/grid/unstructured_mesh.hpp"

#include "pytensor_utils.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


namespace fastscapelib
{

    using py_profile_grid = fs::profile_grid_xt<fs::py_selector>;
    using py_raster_grid = fs::raster_grid_xt<fs::py_selector, fs::raster_connect::queen>;
    using py_unstructured_mesh = fs::unstructured_mesh_xt<fs::py_selector>;

    template <class G>
    void add_grid_static_properties(py::class_<G>& pyg)
    {
        pyg.def_property_readonly_static("is_structured",
                                         [](py::object /*self*/) { return G::is_structured(); })
            .def_property_readonly_static("is_uniform",
                                          [](py::object /*self*/) { return G::is_uniform(); })
            .def_property_readonly_static("max_neighbors",
                                          [](py::object /*self*/) { return G::max_neighbors(); });
    }

    template <class G>
    class py_grid_funcs
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
            if (idx >= g.size())
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

    class py_raster_grid_funcs
    {
    public:
        using size_type = typename py_raster_grid::size_type;
        using indices_type = typename py_raster_grid::neighbors_indices_raster_type;
        using neighbors_type = typename py_raster_grid::neighbors_raster_type;

        py_raster_grid_funcs() = default;

        std::function<indices_type(py_raster_grid&, size_type, size_type)> m_neighbors_indices
            = [this](py_raster_grid& g, size_type row, size_type col)
        {
            check_in_bounds(g, row, col);
            return g.neighbors_indices(row, col);
        };

        std::function<neighbors_type(py_raster_grid&, size_type, size_type)> m_neighbors
            = [this](py_raster_grid& g, size_type row, size_type col)
        {
            check_in_bounds(g, row, col);
            return g.neighbors(row, col);
        };

    private:
        inline void check_in_bounds(const py_raster_grid& g, size_type row, size_type col)
        {
            const auto& shape = g.shape();
            if (row >= shape[0] || col >= shape[1])
            {
                throw std::out_of_range("grid index out of range");
            }
        }
    };

}

#endif  // PYFASTSCAPELIB_GRID_H
