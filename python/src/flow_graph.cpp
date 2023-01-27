#include <memory>

#include "xtensor-python/pytensor.hpp"

#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/stl.h"

#include "grid.hpp"
#include "flow_graph.hpp"

#include "fastscapelib/flow/flow_operator.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/flow_snapshot.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"

namespace py = pybind11;
namespace fs = fastscapelib;


void
add_flow_graph_bindings(py::module& m)
{
    /*
     * Flow graph implementation
     */

    py::class_<fs::py_flow_graph_impl>(m, "FlowGraphImpl")
        .def_property_readonly("receivers", &fs::py_flow_graph_impl::receivers)
        .def_property_readonly("receivers_count", &fs::py_flow_graph_impl::receivers_count)
        .def_property_readonly("receivers_distance", &fs::py_flow_graph_impl::receivers_distance)
        .def_property_readonly("receivers_weight", &fs::py_flow_graph_impl::receivers_weight)
        .def_property_readonly("donors", &fs::py_flow_graph_impl::donors)
        .def_property_readonly("donors_count", &fs::py_flow_graph_impl::donors_count)
        .def_property_readonly("dfs_indices", &fs::py_flow_graph_impl::dfs_indices)
        .def_property_readonly("basins", &fs::py_flow_graph_impl::basins);

    /*
     * Flow operators
     */

    py::class_<fs::flow_operator, std::shared_ptr<fs::flow_operator>>(m, "FlowOperator");

    py::class_<fs::single_flow_router, fs::flow_operator, std::shared_ptr<fs::single_flow_router>>(
        m, "SingleFlowRouter")
        .def(py::init<>());

    py::class_<fs::pflood_sink_resolver,
               fs::flow_operator,
               std::shared_ptr<fs::pflood_sink_resolver>>(m, "PFloodSinkResolver")
        .def(py::init<>());

    py::enum_<fs::mst_method> mst_method(m, "MSTMethod", py::arithmetic());
    mst_method.value("KRUSKAL", fs::mst_method::kruskal).value("BORUVKA", fs::mst_method::boruvka);

    py::enum_<fs::mst_route_method> mst_route_method(m, "MSTRouteMethod", py::arithmetic());
    mst_route_method.value("BASIC", fs::mst_route_method::basic)
        .value("CARVE", fs::mst_route_method::carve);

    py::class_<fs::mst_sink_resolver, fs::flow_operator, std::shared_ptr<fs::mst_sink_resolver>>(
        m, "MSTSinkResolver")
        .def(py::init<fs::mst_method, fs::mst_route_method>(),
             py::arg("basin_method") = fs::mst_method::kruskal,
             py::arg("route_method") = fs::mst_route_method::carve)
        .def_readwrite("basin_method", &fs::mst_sink_resolver::m_basin_method)
        .def_readwrite("route_method", &fs::mst_sink_resolver::m_route_method);

    py::class_<fs::flow_snapshot, fs::flow_operator, std::shared_ptr<fs::flow_snapshot>>(
        m, "FlowSnapshot")
        .def(py::init<std::string, bool, bool>(),
             py::arg("snapshot_name"),
             py::arg("save_graph") = true,
             py::arg("save_elevation") = false);

    /*
     * Flow graph
     */

    py::class_<fs::py_flow_graph> pyfgraph(m, "FlowGraph");

    fs::register_py_flow_graph_init<fs::py_profile_grid>(pyfgraph);
    fs::register_py_flow_graph_init<fs::py_raster_grid>(pyfgraph);

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

    pyfgraph.def("basins", &fs::py_flow_graph::basins);
}
