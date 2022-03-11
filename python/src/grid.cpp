#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"

#include "fastscapelib/grid.hpp"
#include "fastscapelib/profile_grid.hpp"
#include "fastscapelib/raster_grid.hpp"
#include "fastscapelib/unstructured_mesh.hpp"

#include "pytensor_utils.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


void
add_grid_bindings(py::module& m)
{
    // ==== Binding of the node_status enumeration ==== //
    py::enum_<fs::node_status> node_status(
        m,
        "NodeStatus",
        py::arithmetic(),
        "Status of grid/mesh nodes either inside the domain or on the domain boundary.");
    node_status.value("CORE", fs::node_status::core)
        .value("FIXED_VALUE_BOUNDARY", fs::node_status::fixed_value_boundary)
        .value("FIXED_GRADIENT_BOUNDARY", fs::node_status::fixed_gradient_boundary)
        .value("LOOPED_BOUNDARY", fs::node_status::looped_boundary);

    // ==== Binding of the node structure ==== //
    py::class_<fs::node>(m, "Node")
        .def(py::init<std::size_t, fs::node_status>())
        .def_readwrite("idx", &fs::node::idx, "Node index.")
        .def_readwrite("status", &fs::node::status, "Node status.");

    // ==== Binding of the neighbor structure ==== //
    py::class_<fs::neighbor>(m, "Neighbor")
        .def(py::init<std::size_t, double, fs::node_status>())
        .def("__eq__", &fs::neighbor::operator==)
        .def_readwrite("idx", &fs::neighbor::idx, "Neighbor index.")
        .def_readwrite("distance", &fs::neighbor::distance, "Neighbor distance.")
        .def_readwrite("status", &fs::neighbor::status, "Neighbor status.");

    // ==== Binding of the UnstructuredMesh class ==== //
    using unstructured_mesh_py = fs::unstructured_mesh_xt<fs::pyarray_selector>;
    py::class_<unstructured_mesh_py> umesh(m, "UnstructuredMesh");
    umesh.def(py::init());

    // umesh.def_property_readonly("size", [](const unstructured_mesh_py& mesh) { return
    // mesh.size(); })
    //      .def_property_readonly("status_at_nodes", [](const unstructured_mesh_py& mesh) { return
    //      mesh.status_at_nodes(); }) .def("neighbors", [](unstructured_mesh_py& mesh, std::size_t
    //      idx) { return mesh.neighbors(idx); });

    // ==== Binding of the profile_boundary_status class ==== //
    py::class_<fs::profile_boundary_status>(m, "ProfileBoundaryStatus")
        .def(py::init<const fs::node_status>())
        .def(py::init<const fs::node_status, const fs::node_status>())
        .def(py::init<const std::array<fs::node_status, 2>&>())

        .def_readwrite("left", &fs::profile_boundary_status::left, "Left boundary status.")
        .def_readwrite("right", &fs::profile_boundary_status::right, "Right boundary status.")
        .def_property_readonly("is_horizontal_looped",
                               &fs::profile_boundary_status::is_horizontal_looped,
                               "Horizontal looped status.");

    // ==== Binding of the ProfileGrid class ==== //
    py::class_<fs::profile_grid> pgrid(m, "ProfileGrid");
    pgrid.def(py::init<fs::profile_grid::size_type,
                       fs::profile_grid::spacing_type,
                       const fs::profile_grid::boundary_status_type&,
                       const std::vector<fs::node>>());
    pgrid.def(py::init(
        [](std::size_t size,
           fs::profile_grid::spacing_type spacing,
           const std::array<fs::node_status, 2>& bs,
           const std::vector<std::pair<std::size_t, fs::node_status>>& ns)
        {
            std::vector<fs::node> node_vec;
            for (auto&& node : ns)
            {
                node_vec.push_back({ node.first, node.second });
            }
            return std::make_unique<fs::profile_grid>(size, spacing, bs, node_vec);
        }));

    pgrid.def_static("from_length", &fs::profile_grid::from_length);

    pgrid.def_property_readonly("size", [](const fs::profile_grid& g) { return g.size(); })
        .def_property_readonly("shape", [](const fs::profile_grid& g) { return g.shape(); })
        .def_property_readonly("spacing", [](const fs::profile_grid& g) { return g.spacing(); })
        .def_property_readonly("length", [](const fs::profile_grid& g) { return g.length(); })
        .def_property_readonly("status_at_nodes",
                               [](const fs::profile_grid& g) { return g.status_at_nodes(); })
        .def("neighbors", [](fs::profile_grid& g, std::size_t idx) { return g.neighbors(idx); });

    // ==== Binding of the raster_node structure ==== //
    py::class_<fs::raster_node>(m, "RasterNode")
        .def(py::init<std::size_t, std::size_t, fs::node_status>())
        .def_readwrite("row", &fs::raster_node::row, "Node row index.")
        .def_readwrite("col", &fs::raster_node::col, "Node column index.")
        .def_readwrite("status", &fs::raster_node::status, "Node status.");

    // ==== Binding of the raster_neighbor structure ==== //
    py::class_<fs::raster_neighbor>(m, "RasterNeighbor")
        .def(py::init<std::size_t, std::size_t, std::size_t, double, fs::node_status>())
        .def("__eq__", &fs::raster_neighbor::operator==)
        .def_readwrite("flatten_idx", &fs::raster_neighbor::flatten_idx, "Neighbor flatten index.")
        .def_readwrite("row", &fs::raster_neighbor::row, "Neighbor row index.")
        .def_readwrite("col", &fs::raster_neighbor::col, "Neighbor col index.")
        .def_readwrite("distance", &fs::raster_neighbor::distance, "Neighbor distance.")
        .def_readwrite("status", &fs::raster_neighbor::status, "Neighbor status.");

    // ==== Binding of the raster_boundary_status structure ==== //
    py::class_<fs::raster_boundary_status>(m, "RasterBoundaryStatus")
        .def(py::init<const fs::node_status>())
        .def(py::init<const std::array<fs::node_status, 4>&>())

        .def_readwrite("left", &fs::raster_boundary_status::left, "Left boundary status.")
        .def_readwrite("right", &fs::raster_boundary_status::right, "Right boundary status.")
        .def_readwrite("bottom", &fs::raster_boundary_status::bottom, "Bottom boundary status.")
        .def_readwrite("top", &fs::raster_boundary_status::top, "Top boundary status.")
        .def_property_readonly("is_horizontal_looped",
                               &fs::raster_boundary_status::is_horizontal_looped,
                               "Horizontal periodicity.")
        .def_property_readonly("is_vertical_looped",
                               &fs::raster_boundary_status::is_vertical_looped,
                               "Vertical periodicity.");

    // ==== Binding of the RasterGrid class ==== //
    py::class_<fs::raster_grid> rgrid(m, "RasterGrid");
    rgrid.def(py::init<const fs::raster_grid::shape_type&,
                       const xt::pytensor<double, 1>&,
                       const fs::raster_grid::boundary_status_type&,
                       const std::vector<fs::raster_node>>());

    rgrid.def_static("from_length",
                     [](const fs::raster_grid::shape_type& shape,
                        const xt::pytensor<double, 1>& length,
                        const fs::raster_grid::boundary_status_type& boundary,
                        const std::vector<fs::raster_node> status)
                     { return fs::raster_grid::from_length(shape, length, boundary, status); });

    rgrid.def_property_readonly("size", [](const fs::raster_grid& g) { return g.size(); })
        .def_property_readonly("shape", [](const fs::raster_grid& g) { return g.shape(); })
        .def_property_readonly("spacing",
                               [](const fs::raster_grid& g) -> xt::pytensor<double, 1>
                               { return g.spacing(); })
        .def_property_readonly("length",
                               [](const fs::raster_grid& g) -> xt::pytensor<double, 1>
                               { return g.length(); })
        .def_property_readonly("status_at_nodes",
                               [](const fs::raster_grid& g) { return g.status_at_nodes(); });
}
