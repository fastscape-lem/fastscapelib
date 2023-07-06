#ifndef PYFASTSCAPELIB_GRID_H
#define PYFASTSCAPELIB_GRID_H

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "pybind11/pybind11.h"

#include "xtensor-python/pytensor.hpp"

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/grid/profile_grid.hpp"
#include "fastscapelib/grid/raster_grid.hpp"
#include "fastscapelib/grid/trimesh.hpp"

#include "pytensor_utils.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


namespace fastscapelib
{

    using py_profile_grid = fs::profile_grid_xt<fs::py_selector>;
    using py_raster_grid = fs::raster_grid_xt<fs::py_selector, fs::raster_connect::queen>;
    using py_trimesh = fs::trimesh_xt<fs::py_selector>;

    template <class G>
    void register_grid_static_properties(py::class_<G>& pyg)
    {
        pyg.def_property_readonly_static(
               "is_structured",
               [](py::object /*self*/) { return G::is_structured(); },
               "Returns True if the grid type is structured.")
            .def_property_readonly_static(
                "is_uniform",
                [](py::object /*self*/) { return G::is_uniform(); },
                "Returns True if the grid type has uniform node spacing.")
            .def_property_readonly_static(
                "n_neighbors_max",
                [](py::object /*self*/) { return G::n_neighbors_max(); },
                "Returns the maximum number of grid node neighbors.");
    }

    /*
    ** Use lambdas since directly using method pointers doesn't seem to
    ** always work (issues with GCC 11)
    */
    template <class G>
    void register_base_grid_properties(py::class_<G>& pyg)
    {
        pyg.def_property_readonly(
            "size", [](const G& g) { return g.size(); }, "Returns the total number of grid nodes.");
        pyg.def_property_readonly(
            "shape",
            [](const G& g) { return g.shape(); },
            "Returns the shape of the grid node arrays.");
    }

    /*
    ** Use lambdas since directly using method pointers doesn't seem to
    ** always work (issues with GCC 11)
    */
    template <class G>
    void register_structured_grid_properties(py::class_<G>& pyg)
    {
        static_assert(G::is_structured(),
                      "cannot register structured grid properties to a non-structured grid type");

        using property_t = std::
            conditional_t<std::is_same<G, py_raster_grid>::value, xt::pytensor<double, 1>, double>;

        pyg.def_property_readonly(
            "spacing",
            [](const G& g) -> property_t { return g.spacing(); },
            "Returns the fixed distance between two adjacent nodes for each dimension.");
        pyg.def_property_readonly(
            "length",
            [](const G& g) -> property_t { return g.length(); },
            "Returns the total length of the grid for each of its dimension.");
    }

    template <class G>
    class py_grid_funcs
    {
    public:
        using size_type = typename G::size_type;
        using count_type = typename G::size_type;
        using indices_type = typename G::neighbors_indices_type;
        using distances_type = typename G::neighbors_distances_type;
        using neighbors_type = typename G::neighbors_type;

        py_grid_funcs() = default;

        std::function<xt::pytensor<size_type, 1>(const G&)> m_nodes_indices = [this](const G& g)
        {
            auto indices = g.nodes_indices();
            auto output = xt::pyarray<size_type>::from_shape({ g.size() });
            std::copy(indices.begin(), indices.end(), output.begin());

            return std::move(output);
        };

        std::function<xt::pytensor<size_type, 1>(const G&, node_status)> m_nodes_indices2
            = [this](const G& g, node_status status)
        {
            auto indices = g.nodes_indices(status);

            // need a temporary vector as the container size is not known in advance
            // and pytensor doesn't support well resizing
            std::vector<size_type> temp;

            for (auto& idx : g.nodes_indices(status))
            {
                temp.push_back(idx);
            }

            auto output = xt::pyarray<size_type>::from_shape({ temp.size() });
            std::copy(temp.begin(), temp.end(), output.begin());

            return std::move(output);
        };

        std::function<xt::pyarray<node_status>(const G&)> m_nodes_status
            = [this](const G& g) { return g.nodes_status(); };

        std::function<node_status(const G&, size_type)> m_nodes_status2
            = [this](const G& g, size_type idx)
        {
            check_in_bounds(g, idx);
            return g.nodes_status(idx);
        };

        std::function<double(const G&, size_type)> m_nodes_areas = [this](const G& g, size_type idx)
        {
            check_in_bounds(g, idx);
            return g.nodes_areas(idx);
        };

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
    void register_grid_methods(py::class_<G>& pyg)
    {
        auto grid_funcs = py_grid_funcs<G>();

        pyg.def("nodes_indices",
                grid_funcs.m_nodes_indices,
                R"doc(nodes_indices(*args) -> numpy.ndarray

            Returns the (flattened) indices of grid nodes, possibly filtered
            by node status.

            Overloaded method that supports the following signatures:

            1. ``nodes_indices() -> numpy.ndarray``

            2. ``nodes_indices(status: NodeStatus) -> numpy.ndarray``

            Parameters
            ----------
            status : :class:`NodeStatus`
                The status of the grid nodes to filter in.

            Returns
            -------
            indices : numpy.ndarray
                A new, 1-dimensional array with grid node (filtered) flat indices.

            )doc");

        pyg.def("nodes_indices", grid_funcs.m_nodes_indices2, py::arg("status"));

        pyg.def("nodes_status",
                grid_funcs.m_nodes_status,
                R"doc(nodes_status(*args) -> NodeStatus | numpy.ndarray

            Returns the status of a given node or all nodes.

            Overloaded method that supports the following signatures:

            1. ``nodes_status() -> numpy.ndarray``

            2. ``nodes_status(idx: int) -> NodeStatus``

            Parameters
            ----------
            idx : int
                The grid node (flat) index.

            Returns
            -------
            status : :class:`NodeStatus` or numpy.ndarray
                If an array is returned, it is a copy of the original array
                and has the :class:`NodeStatus` enum class values (int) as
                values.

            )doc");

        pyg.def("nodes_status", grid_funcs.m_nodes_status2, py::arg("idx"));

        pyg.def("nodes_areas",
                py::overload_cast<>(&G::nodes_areas, py::const_),
                R"doc(nodes_areas(*args) -> float | numpy.ndarray

            Returns the area of the direct vicinity of a given node or all
            nodes.

            Overloaded method that supports the following signatures:

            1. ``nodes_areas() -> numpy.ndarray``

            2. ``nodes_areas(idx: int) -> float``

            Parameters
            ----------
            idx : int
                The grid node (flat) index.

            )doc");

        pyg.def("nodes_areas", grid_funcs.m_nodes_areas, py::arg("idx"));

        pyg.def("neighbors_count",
                grid_funcs.m_neighbors_count,
                py::arg("idx"),
                R"doc(neighbors_count(idx: int) -> int

            Returns the number of neighbors at a given node.

            Parameters
            ----------
            idx : int
                The grid node (flat) index.

            )doc");

        pyg.def("neighbors_indices",
                grid_funcs.m_neighbors_indices,
                py::arg("idx"),
                R"doc(neighbors_indices(idx: int) -> numpy.ndarray

            Returns the (flat) indices of the neighbors of a given node.

            Parameters
            ----------
            idx : int
                The grid node (flat) index.

            Returns
            -------
            indices : numpy.ndarray
                A 1-dimensional array of node (flat) indices (int).

            Notes
            -----
            For a :class:`RasterGrid`, it is also possible to provide
            ``row`` and ``col`` indices as arguments.

            )doc");

        pyg.def("neighbors_distances",
                grid_funcs.m_neighbors_distances,
                py::arg("idx"),
                R"doc(neighbors_distances(idx: int) -> numpy.ndarray

            Returns the distances from a given node to its neighbors.

            Parameters
            ----------
            idx : int
                The grid node (flat) index.

            Returns
            -------
            distances : numpy.ndarray
                A 1-dimensional array of node distances (float). The neighbor order is
                the same than the one returned by ``neighbors_indices(idx)``.

            )doc");

        pyg.def("neighbors",
                grid_funcs.m_neighbors,
                py::arg("idx"),
                R"doc(neighbors(idx: int) -> List[Neighbor]

            Returns the neighbors of a given node.

            Parameters
            ----------
            idx : int
                The grid node (flat) index.

            Returns
            -------
            neighbors : list
                A list of :class:`Neighbor` instances.

            Notes
            -----
            For a :class:`RasterGrid`, it is also possible to provide
            ``row`` and ``col`` indices as arguments, which returns a
            list of :class:`RasterNeighbor` instances.

            )doc");
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
