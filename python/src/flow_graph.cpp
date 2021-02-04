#include "flow_graph.hpp"
#include "flow_router.hpp"
#include "sink_resolver.hpp"

#include "pybind11/pybind11.h"


namespace py = pybind11;
namespace fs = fastscapelib;


void add_flow_graph_bindings(py::module& m) {

    using namespace fs::detail;

    // ==== Binding of the FlowGraph class ==== //
    py::class_<flow_graph_facade>(m, "FlowGraph")
        .def(py::init<fs::profile_grid&, std::shared_ptr<flow_router_method>, std::shared_ptr<sink_resolver_method> >(), 
                        py::arg("grid"), py::arg("flow_router"), py::arg("sink_resolver") = std::make_shared<no_sink_resolver_method>())
        .def(py::init<fs::raster_grid&, std::shared_ptr<flow_router_method>, std::shared_ptr<sink_resolver_method> >(), 
                        py::arg("grid"), py::arg("flow_router"), py::arg("sink_resolver") = std::make_shared<no_sink_resolver_method>())
        .def("update_routes", &flow_graph_facade::update_routes)
        .def("receivers", &flow_graph_facade::receivers)
        .def("receivers_count", &flow_graph_facade::receivers_count)
        .def("receivers_distance", &flow_graph_facade::receivers_distance)
        .def("receivers_weight", &flow_graph_facade::receivers_weight)
        .def("donors", &flow_graph_facade::donors)
        .def("donors_count", &flow_graph_facade::donors_count)
        .def("dfs_stack", &flow_graph_facade::dfs_stack)
        .def("accumulate", py::overload_cast<const xt::pyarray<double>&>(&flow_graph_facade::accumulate, py::const_))
        .def("accumulate", py::overload_cast<const double&>(&flow_graph_facade::accumulate, py::const_));
}
 