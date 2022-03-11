#include <memory>

#include "pybind11/pybind11.h"

#include "fastscapelib/profile_grid.hpp"
#include "fastscapelib/raster_grid.hpp"

#include "flow_graph.hpp"
#include "flow_router.hpp"
#include "sink_resolver.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


void
add_flow_graph_bindings(py::module& m)
{
    using flow_graph = fs::detail::flow_graph_facade;

    // ==== Binding of the FlowGraph class ==== //
    py::class_<flow_graph>(m, "FlowGraph")
        .def(py::init<fs::profile_grid&,
                      std::shared_ptr<fs::detail::flow_router_method>,
                      std::shared_ptr<fs::detail::sink_resolver_method>>(),
             py::arg("grid"),
             py::arg("flow_router"),
             py::arg("sink_resolver") = std::make_shared<fs::detail::no_sink_resolver_method>())
        .def(py::init<fs::raster_grid&,
                      std::shared_ptr<fs::detail::flow_router_method>,
                      std::shared_ptr<fs::detail::sink_resolver_method>>(),
             py::arg("grid"),
             py::arg("flow_router"),
             py::arg("sink_resolver") = std::make_shared<fs::detail::no_sink_resolver_method>())
        .def("update_routes", &flow_graph::update_routes)
        .def("receivers", &flow_graph::receivers)
        .def("receivers_count", &flow_graph::receivers_count)
        .def("receivers_distance", &flow_graph::receivers_distance)
        .def("receivers_weight", &flow_graph::receivers_weight)
        .def("donors", &flow_graph::donors)
        .def("donors_count", &flow_graph::donors_count)
        .def("dfs_stack", &flow_graph::dfs_stack)
        .def("accumulate",
             py::overload_cast<const xt::pyarray<double, xt::layout_type::row_major>&>(
                 &flow_graph::accumulate, py::const_))
        .def("accumulate", py::overload_cast<const double&>(&flow_graph::accumulate, py::const_));
}
