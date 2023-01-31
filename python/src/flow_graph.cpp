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
        .def_property_readonly("single_flow", &fs::py_flow_graph_impl::single_flow)
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

    py::enum_<fs::flow_direction> flow_direction(
        m, "FlowDirection", py::arithmetic(), "Flow direction strategy.");
    flow_direction.value("UNDEFINED", fs::flow_direction::undefined)
        .value("SINGLE", fs::flow_direction::single)
        .value("MULTIPLE", fs::flow_direction::multiple);

    py::class_<fs::flow_operator, std::shared_ptr<fs::flow_operator>>(m, "FlowOperator")
        .def_property_readonly("name", &fs::flow_operator::name)
        .def("_repr_inline_", [](const fs::flow_operator&) { return ""; })
        .def("__repr__",
             [](const fs::flow_operator& op)
             {
                 auto class_name = py::cast(op).get_type().attr("__name__").cast<std::string>();
                 return "class_name (" + op.name() + ")";
             });

    py::class_<fs::single_flow_router, fs::flow_operator, std::shared_ptr<fs::single_flow_router>>
        srouter_op(m, "SingleFlowRouter");

    srouter_op.def(py::init<>());
    fs::register_operator_static_attrs(srouter_op);
    srouter_op.def("__repr__", [](const fs::single_flow_router& op) { return "SingleFlowRouter"; });

    py::class_<fs::pflood_sink_resolver,
               fs::flow_operator,
               std::shared_ptr<fs::pflood_sink_resolver>>
        pflood_op(m, "PFloodSinkResolver");

    pflood_op.def(py::init<>());
    fs::register_operator_static_attrs(pflood_op);
    pflood_op.def("__repr__",
                  [](const fs::pflood_sink_resolver& op) { return "PFloodSinkResolver"; });

    py::enum_<fs::mst_method> mst_method(m, "MSTMethod", py::arithmetic());
    mst_method.value("KRUSKAL", fs::mst_method::kruskal).value("BORUVKA", fs::mst_method::boruvka);

    py::enum_<fs::mst_route_method> mst_route_method(m, "MSTRouteMethod", py::arithmetic());
    mst_route_method.value("BASIC", fs::mst_route_method::basic)
        .value("CARVE", fs::mst_route_method::carve);

    py::class_<fs::mst_sink_resolver, fs::flow_operator, std::shared_ptr<fs::mst_sink_resolver>>
        mst_op(m, "MSTSinkResolver");

    mst_op.def(py::init<fs::mst_method, fs::mst_route_method>(),
               py::arg("basin_method") = fs::mst_method::kruskal,
               py::arg("route_method") = fs::mst_route_method::carve);
    fs::register_operator_static_attrs(mst_op);
    mst_op.def_readwrite("basin_method", &fs::mst_sink_resolver::m_basin_method);
    mst_op.def_readwrite("route_method", &fs::mst_sink_resolver::m_route_method);
    mst_op.def("__repr__",
               [](const fs::mst_sink_resolver& op)
               {
                   std::string bmeth
                       = op.m_basin_method == fs::mst_method::kruskal ? "kruskal" : "boruvka";
                   std::string rmeth
                       = op.m_route_method == fs::mst_route_method::basic ? "basic" : "carve";
                   return "MSTSinkResolver (basin=" + bmeth + ", route=" + rmeth + ")";
               });

    py::class_<fs::flow_snapshot, fs::flow_operator, std::shared_ptr<fs::flow_snapshot>>
        snapshot_op(m, "FlowSnapshot");
    fs::register_operator_static_attrs(snapshot_op);
    snapshot_op.def(py::init<std::string, bool, bool>(),
                    py::arg("snapshot_name"),
                    py::arg("save_graph") = true,
                    py::arg("save_elevation") = false);
    snapshot_op.def_property_readonly("snapshot_name", &fs::flow_snapshot::snapshot_name);
    snapshot_op.def_property_readonly("save_graph", &fs::flow_snapshot::save_graph);
    snapshot_op.def_property_readonly("save_elevation", &fs::flow_snapshot::save_elevation);
    snapshot_op.def("__repr__",
                    [](const fs::flow_snapshot& op)
                    {
                        std::string save_graph = op.save_graph() == true ? "True" : "False";
                        std::string save_elevation = op.save_elevation() == true ? "True" : "False";
                        auto repr = "FlowSnapshot '" + op.snapshot_name() + "' ";
                        repr += "(graph=" + save_graph + ", ";
                        repr += "elevation=" + save_elevation + ")";
                        return repr;
                    });

    /*
     * Flow graph
     */

    py::class_<fs::py_flow_graph> pyfgraph(m, "FlowGraph");

    fs::register_py_flow_graph_init<fs::py_profile_grid>(pyfgraph);
    fs::register_py_flow_graph_init<fs::py_raster_grid>(pyfgraph);
    fs::register_py_flow_graph_init<fs::py_unstructured_mesh>(pyfgraph);

    pyfgraph.def_property_readonly("operators", &fs::py_flow_graph::operators);
    pyfgraph.def("impl", &fs::py_flow_graph::impl, py::return_value_policy::reference);
    pyfgraph.def_property_readonly("graph_snapshot_keys", &fs::py_flow_graph::graph_snapshot_keys);
    pyfgraph.def("graph_snapshot", &fs::py_flow_graph::graph_snapshot);
    pyfgraph.def_property_readonly("elevation_snapshot_keys",
                                   &fs::py_flow_graph::elevation_snapshot_keys);
    pyfgraph.def("elevation_snapshot", &fs::py_flow_graph::elevation_snapshot);
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

    pyfgraph.def("__repr__",
                 [](const fs::py_flow_graph& pyfg)
                 {
                     auto nnodes = std::to_string(pyfg.size());
                     std::string repr = "<FlowGraph (" + nnodes + " nodes)>";
                     repr += "\nOperators:";
                     for (auto& op : pyfg.operators())
                     {
                         repr += "\n    " + py::repr(op).cast<std::string>();
                     }
                     if (pyfg.operators().empty())
                     {
                         repr += "\n    *empty*";
                     }
                     repr += "\n";

                     return repr;
                 });
}
