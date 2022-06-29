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
    py::class_<fs::py_flow_graph_impl>(m, "FlowGraphImpl")
        .def("receivers", &fs::py_flow_graph_impl::receivers)
        .def("receivers_count", &fs::py_flow_graph_impl::receivers_count)
        .def("receivers_distance", &fs::py_flow_graph_impl::receivers_distance)
        .def("receivers_weight", &fs::py_flow_graph_impl::receivers_weight)
        .def("donors", &fs::py_flow_graph_impl::donors)
        .def("donors_count", &fs::py_flow_graph_impl::donors_count)
        .def("dfs_stack", &fs::py_flow_graph_impl::dfs_stack);

    py::class_<fs::py_flow_graph> pyfgraph(m, "FlowGraph");

    fs::add_init_methods<fs::py_profile_grid>(pyfgraph);
    fs::add_init_methods<fs::py_raster_grid>(pyfgraph);

    pyfgraph.def("impl", &fs::py_flow_graph::impl, py::return_value_policy::reference);

    pyfgraph.def("update_routes", &fs::py_flow_graph::update_routes);

    using data_array_type = fs::py_flow_graph::data_array_type;
    using data_type = fs::py_flow_graph::data_type;

    pyfgraph
        .def("accumulate",
             py::overload_cast<data_array_type&, const data_array_type&>(
                 &fs::py_flow_graph::accumulate, py::const_))
        .def("accumulate",
             py::overload_cast<data_array_type&, data_type>(&fs::py_flow_graph::accumulate,
                                                            py::const_))
        .def("accumulate",
             py::overload_cast<const data_array_type&>(&fs::py_flow_graph::accumulate, py::const_))
        .def("accumulate",
             py::overload_cast<data_type>(&fs::py_flow_graph::accumulate, py::const_));
}
