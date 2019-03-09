#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "fastscapelib/utils.hpp"
#include "fastscapelib/grid.hpp"

#include "xtensor_utils.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


void add_grid_bindings(py::module& m);


void add_grid_bindings(py::module& m) {

    using profile_grid_py = fs::profile_grid_xt<fs::pytensor_selector>;

    py::enum_<fs::node_status>(m, "NodeStatus", py::arithmetic(),
                               "Status of grid/mesh nodes either inside the domain "
                               "or on the domain boundary.")
        .value("CORE", fs::node_status::core)
        .value("FIXED_VALUE_BOUNDARY", fs::node_status::fixed_value_boundary)
        .value("FIXED_GRADIENT_BOUNDARY", fs::node_status::fixed_gradient_boundary)
        .value("LOOPED_BOUNDARY", fs::node_status::looped_boundary);

    py::class_<fs::edge_nodes_status>(m, "EdgeNodesStatus")
        .def(py::init<fs::node_status, fs::node_status>())
        .def_readonly("left", &fs::edge_nodes_status::left)
        .def_readonly("right", &fs::edge_nodes_status::right);

    py::class_<fs::node>(m, "Node")
        .def(py::init<std::size_t, fs::node_status>())
        .def_readonly("idx", &fs::node::idx)
        .def_readonly("status", &fs::node::status);

    py::class_<profile_grid_py>(m, "ProfileGrid")
        .def(py::init<std::size_t, double, const fs::edge_nodes_status,
                      const std::vector<fs::node>&>())
        .def_property_readonly("size", &profile_grid_py::size)
        .def_property_readonly("spacing", &profile_grid_py::spacing)
        .def_property_readonly("node_status", &profile_grid_py::node_status);

}
