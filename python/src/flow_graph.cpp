#include <memory>

#include "xtensor-python/pytensor.hpp"

#include "pybind11/pybind11.h"

#include "grid.hpp"
#include "flow_graph.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


void
add_flow_graph_bindings(py::module& m)
{
    using py_flow_graph = fs::detail::py_flow_graph;

    py::class_<py_flow_graph>(m, "FlowGraph")
        .def(py::init(
                 [](fs::py_profile_grid& grid,
                    const std::shared_ptr<fs::flow_router> ptr_router,
                    const std::shared_ptr<fs::sink_resolver> ptr_resolver)
                 {
                     auto router = std::static_pointer_cast<fs::single_flow_router>(ptr_router);
                     return std::make_unique<py_flow_graph>(grid, *router, *ptr_resolver);
                 }),
             py::arg("grid"),
             py::arg("flow_router"),
             py::arg("sink_resolver"))
        .def(py::init(
                 [](fs::py_raster_grid& grid,
                    const std::shared_ptr<fs::flow_router> ptr_router,
                    const std::shared_ptr<fs::sink_resolver> ptr_resolver)
                 {
                     auto router = dynamic_cast<fs::single_flow_router&>(*ptr_router);
                     return std::make_unique<py_flow_graph>(grid, router, *ptr_resolver);
                 }),
             py::arg("grid"),
             py::arg("flow_router"),
             py::arg("sink_resolver"))
        .def("update_routes", &py_flow_graph::update_routes)
        .def("receivers", &py_flow_graph::receivers)
        .def("receivers_count", &py_flow_graph::receivers_count)
        .def("receivers_distance", &py_flow_graph::receivers_distance)
        .def("receivers_weight", &py_flow_graph::receivers_weight)
        .def("donors", &py_flow_graph::donors)
        .def("donors_count", &py_flow_graph::donors_count)
        .def("dfs_stack", &py_flow_graph::dfs_stack)
        .def("accumulate",
             py::overload_cast<const xt::pyarray<double, xt::layout_type::row_major>&>(
                 &py_flow_graph::accumulate, py::const_))
        .def("accumulate",
             py::overload_cast<const double&>(&py_flow_graph::accumulate, py::const_));
}
