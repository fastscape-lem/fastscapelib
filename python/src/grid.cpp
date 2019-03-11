#include <array>
#include <memory>
#include <utility>
#include <vector>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "fastscapelib/utils.hpp"
#include "fastscapelib/grid.hpp"

#include "xtensor_utils.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


using profile_grid_py = fs::profile_grid_xt<fs::pytensor_selector>;

void add_grid_bindings(py::module& m);


void add_grid_bindings(py::module& m) {

    py::enum_<fs::node_status>(m, "NodeStatus", py::arithmetic(),
                               "Status of grid/mesh nodes either inside the domain "
                               "or on the domain boundary.")
        .value("CORE", fs::node_status::core)
        .value("FIXED_VALUE_BOUNDARY", fs::node_status::fixed_value_boundary)
        .value("FIXED_GRADIENT_BOUNDARY", fs::node_status::fixed_gradient_boundary)
        .value("LOOPED_BOUNDARY", fs::node_status::looped_boundary);

    py::class_<profile_grid_py>(m, "ProfileGrid")
        .def(py::init(
                 [](std::size_t size,
                    double spacing,
                    const std::array<fs::node_status, 2>& es,
                    const std::vector<std::pair<std::size_t, fs::node_status>>& ns)
                 {
                     std::vector<fs::node> vec;
                     for (auto&& p : ns)
                     {
                         vec.push_back({std::get<0>(p), std::get<1>(p)});
                     }
                     return std::make_unique<profile_grid_py>(size, spacing, es, vec);
                 }))
        .def_property_readonly("size", &profile_grid_py::size)
        .def_property_readonly("spacing", &profile_grid_py::spacing)
        .def_property_readonly("status_at_nodes", &profile_grid_py::status_at_nodes);
}
