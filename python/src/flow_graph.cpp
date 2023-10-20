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
    py::options options;
    options.disable_function_signatures();

    /*
     * Flow graph implementation
     */

    py::class_<fs::py_flow_graph_impl> flow_graph_impl(
        m,
        "FlowGraphImpl",
        R"doc(Flow graph fixed-shape array implementation.

        This class has no initializer. Instead, use the :meth:`FlowGraph.impl` method
        to access an implementation instance from a flow graph instance.

        This class exposes low-level arrays of fixed-shape used to represent the flow
        graph. In general, :class:`FlowGraph` API should be used instead when possible.

        IMPORTANT NOTES:

        This class is not considered as part of Fastscapelib's public API!

        The (read-only) properties of this class return the original arrays (no
        copy). Directly changing value items may invalidate the flow graph.

        )doc");

    flow_graph_impl
        .def_property_readonly("single_flow",
                               &fs::py_flow_graph_impl::single_flow,
                               "Returns True if the graph has single flow directions (flow tree).")
        .def_property_readonly("receivers",
                               &fs::py_flow_graph_impl::receivers,
                               R"doc(Returns the array of flow receiver grid node (flat) indices.

            The shape of the array is either ``(size, 1)`` for single flow graphs or
            ``(size, n_neighbors_max)`` for multiple flow graphs.

            For base-level or pit nodes: flow receiver index = node index.

            Note: values returned by ``receivers[i, j]`` where ``j >= receivers_count[i]``
            are meaningless.

            )doc")
        .def_property_readonly(
            "receivers_count",
            &fs::py_flow_graph_impl::receivers_count,
            "Returns the total number of receivers for each node (1-dimensional array).")
        .def_property_readonly(
            "receivers_distance",
            &fs::py_flow_graph_impl::receivers_distance,
            R"doc(Returns the array of distances from a node to each of its receivers.

            The shape of the array is either ``(size, 1)`` for single flow graphs or
            ``(size, n_neighbors_max)`` for multiple flow graphs.

            Note: values returned by ``receivers_distance[i, j]`` where ``j >= receivers_count[i]``
            are meaningless.

            )doc")
        .def_property_readonly(
            "receivers_weight",
            &fs::py_flow_graph_impl::receivers_weight,
            R"doc(Returns the array of flow partition weights from a node to each of its receivers.

            The shape of the array is either ``(size, 1)`` for single flow graphs or
            ``(size, n_neighbors_max)`` for multiple flow graphs.

            Note: values returned by ``receivers_weight[i, j]`` where ``j >= receivers_count[i]``
            are meaningless.

            )doc")
        .def_property_readonly("donors",
                               &fs::py_flow_graph_impl::donors,
                               R"doc(Returns the array of flow donor grid node (flat) indices.

            The array has the shape ``(size, n_neighbors_max)``.

            Note: values returned by ``donors[i, j]`` where ``j >= donors_count[i]``
            are meaningless.

            )doc")
        .def_property_readonly(
            "donors_count",
            &fs::py_flow_graph_impl::donors_count,
            "Returns the total number of donors for each node (1-dimensional array).")
        .def_property_readonly("dfs_indices",
                               &fs::py_flow_graph_impl::dfs_indices,
                               "Deprecated alias of :attr:`FlowGraphImpl.nodes_indices_bottomup`.")
        .def_property_readonly(
            "nodes_indices_bottomup",
            &fs::py_flow_graph_impl::dfs_indices,
            R"doc(Returns the node indices ordered topologically from base level nodes up to top nodes.

            The order results from graph traversal using depth-first search.
            )doc")
        .def_property_readonly(
            "basins",
            &fs::py_flow_graph_impl::basins,
            R"doc(Returns an array with the id of the basin (catchment) to which each node belongs.

            Array values are valid only for a single-flow graph!

            )doc");

    /*
     * Flow operators
     */

    py::enum_<fs::flow_direction> flow_direction(
        m, "FlowDirection", py::arithmetic(), "Flow direction strategy.");
    flow_direction
        .value(
            "UNDEFINED", fs::flow_direction::undefined, "Pass-through or flow direction agnostic.")
        .value("SINGLE", fs::flow_direction::single, "Single flow (one receiver per node).")
        .value(
            "MULTI", fs::flow_direction::multi, "Multiple flow (one or more receivers per node).");

    py::class_<fs::flow_operator, std::shared_ptr<fs::flow_operator>> flow_op(
        m,
        "FlowOperator",
        R"doc(Flow operator base (abstract) class.

        This class has no initializer. Use insteads its derived classes.

        )doc");

    flow_op.def_property_readonly("name", &fs::flow_operator::name, "Flow operator's name.");

    py::class_<fs::single_flow_router, fs::flow_operator, std::shared_ptr<fs::single_flow_router>>
        srouter_op(m,
                   "SingleFlowRouter",
                   R"doc(Single direction flow router operator.

            This flow operator routes all the flow passing through a grid node
            towards its neighbors node of steepest slope.

            On a raster grid with 8-node connectivity, this is equivalent to the
            so-called "D8" algorithm :cite:p:`OCallaghan1984`.

            )doc");

    fs::register_operator_static_attrs(srouter_op);

    srouter_op.def(py::init<>(),
                   R"doc(__init__(self) -> None

        SingleFlowRouter initializer.

        )doc");
    srouter_op.def("__repr__", [](const fs::single_flow_router& op) { return "SingleFlowRouter"; });

    py::class_<fs::multi_flow_router, fs::flow_operator, std::shared_ptr<fs::multi_flow_router>>
        mrouter_op(m,
                   "MultiFlowRouter",
                   R"doc(Multiple direction flow router operator.

            This flow operator partitions the flow passing through a grid node among
            its downslope neighbor nodes. Flow partitioning is proportional to the
            local slope between a node and its neighbors (power relationship with a
            fixed exponent parameter).

            The fraction :math:`f_{i,j}` of flow routed from node :math:`i` to its
            neighbor :math:`j` is given by

            .. math::

               f_{i,j} = \frac{\max (0, S_{i, j}^p)}{\sum_{k \in N} \max (0, S_{i, k}^p)}

            where :math:`p` is the slope exponent parameter, :math:`S_{i, j}` is
            the slope between :math:`i`, :math:`j` and :math:`N` is the set of
            all neighbors of :math:`i`.

            Depending on the value of the slope exponent parameter, this is
            equivalent to the methods described in :cite:t:`Quinn1991` or
            :cite:t:`Holmgren1994`.

            )doc");

    fs::register_operator_static_attrs(mrouter_op);

    mrouter_op.def(py::init<double>(),
                   py::arg("slope_exp") = 1.0,
                   R"doc(__init__(self, slope_exp: float = 1.0) -> None

        MultiFlowRouter initializer.

        Parameters
        ----------
        slope_exp : float
            Flow partition slope exponent.

        )doc");
    mrouter_op.def_readwrite(
        "slope_exp", &fs::multi_flow_router::m_slope_exp, "Flow partition slope exponent.");
    mrouter_op.def("__repr__",
                   [](const fs::multi_flow_router& op) {
                       return "MultiFlowRouter (slope_exp=" + std::to_string(op.m_slope_exp) + ")";
                   });

    py::class_<fs::pflood_sink_resolver,
               fs::flow_operator,
               std::shared_ptr<fs::pflood_sink_resolver>>
        pflood_op(m,
                  "PFloodSinkResolver",
                  R"doc(Priority-flood sink resolver operator.

            This flow operator fills the closed depressions in the topographic
            surface using the priority flood algorithm +epsilon variant :cite:p:`Barnes2014`.
            This variant prevents flat surfaces and hence ensure
            that the flow can be routed towards the outlets without disruption.

            )doc");

    fs::register_operator_static_attrs(pflood_op);

    pflood_op.def(py::init<>(),
                  R"doc(__init__(self) -> None

        PFloodSinkResolver initializer.

        )doc");
    pflood_op.def("__repr__",
                  [](const fs::pflood_sink_resolver& op) { return "PFloodSinkResolver"; });

    py::enum_<fs::mst_method> mst_method(
        m,
        "MSTMethod",
        py::arithmetic(),
        "Basin simplification (minimum spanning tree) algorithm used by :class:`MSTSinkResolver`");
    mst_method
        .value("KRUSKAL",
               fs::mst_method::kruskal,
               "Kruskal algorithm (simple implementation, O(n log n) complexity)")
        .value("BORUVKA",
               fs::mst_method::boruvka,
               "Boruvka algorithm (complex implementation, O(n) complexity)");

    py::enum_<fs::mst_route_method> mst_route_method(
        m,
        "MSTRouteMethod",
        py::arithmetic(),
        R"doc(Method used by :class:`MSTSinkResolver` for routing the flow
        within inner basins.

        The ``BASIC`` method is the most efficient one but does not result in a
        realistic, planar flow graph.

        The ``CARVE`` method Mimics carving a narrow canyon within the
        depression.

        )doc");

    mst_route_method
        .value("BASIC",
               fs::mst_route_method::basic,
               "Connect the the pit node directly to the depression spill node")
        .value("CARVE",
               fs::mst_route_method::carve,
               "Revert the (unique) flow path between the spill and pit nodes");

    py::class_<fs::mst_sink_resolver, fs::flow_operator, std::shared_ptr<fs::mst_sink_resolver>>
        mst_op(m,
               "MSTSinkResolver",
               R"doc(Minimum Spanning Tree (MST) sink resolver operator.

            This flow operator re-routes the flow trapped in closed depressions
            towards their spill, using an efficient algorithm that explicitly
            computes a graph of (inner and outer) basins and reduces it as a tree
            :cite:p:`Cordonnier2019`.

            It requires a single flow graph as input.

            This operator also use the updated routes in closed depressions to
            fill these with nearly flat surfaces (a tiny slope ensure natural
            flow routing for the operators applied after this one).

            )doc");

    fs::register_operator_static_attrs(mst_op);

    mst_op.def(
        py::init<fs::mst_method, fs::mst_route_method>(),
        py::arg("basin_method") = fs::mst_method::kruskal,
        py::arg("route_method") = fs::mst_route_method::carve,
        R"doc(__init__(self, basin_method: MSTMethod = MSTMethod.KRUSKAL, route_method: MSTRouteMethod = MSTRouteMethod.CARVE) -> None

        MSTSinkResolver initializer.

        Parameters
        ----------
        basin_method : :class:`MSTMethod`
            Basin simplification algorithm selector.
        route_method : :class:`MSTRouteMethod`
            Method used to route the flow within inner basins.

        )doc");
    mst_op.def_readwrite(
        "basin_method", &fs::mst_sink_resolver::m_basin_method, "Basin simplification algorithm.");
    mst_op.def_readwrite("route_method",
                         &fs::mst_sink_resolver::m_route_method,
                         "Method used to route the flow within inner basins.");
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
        snapshot_op(m,
                    "FlowSnapshot",
                    R"doc(Flow snapshot operator.

            A special flow operator used to save intermediate states
            of the flow graph and/or topographic elevation values while
            applying the other operators in chain.

            Those saved states are accessible from
            :meth:`FlowGraph.graph_snapshot` and
            :meth:`FlowGraph.elevation_snapshot`.

            )doc");

    fs::register_operator_static_attrs(snapshot_op);

    snapshot_op.def(
        py::init<std::string, bool, bool>(),
        py::arg("snapshot_name"),
        py::arg("save_graph") = true,
        py::arg("save_elevation") = false,
        R"doc(__init__(self, snapshot_name: str, save_graph: bool = True, save_elevation: bool = False) -> None

        FlowSnapshot initializer.

        Parameters
        ----------
        snapshot_name : str
            Name of the snapshot.
        save_graph : bool
            If True (default), save the intermediate graph state.
        save_elevation : bool
            If True, save the intermediate topography elevation (default: False).

        )doc");
    snapshot_op.def_property_readonly(
        "snapshot_name", &fs::flow_snapshot::snapshot_name, "The snapshot name.");
    snapshot_op.def_property_readonly("save_graph",
                                      &fs::flow_snapshot::save_graph,
                                      "If True, saves the intermediate graph state.");
    snapshot_op.def_property_readonly("save_elevation",
                                      &fs::flow_snapshot::save_elevation,
                                      "If True, saves the intermediate topography elevation.");
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

    py::class_<fs::py_flow_graph> pyfgraph(
        m,
        "FlowGraph",
        "Main class used to compute or follow flow routes on the topographic surface.");

    fs::register_py_flow_graph_init<fs::py_profile_grid>(pyfgraph, true);
    fs::register_py_flow_graph_init<fs::py_raster_grid>(pyfgraph, false);
    fs::register_py_flow_graph_init<fs::py_trimesh>(pyfgraph, false);

    pyfgraph.def_property_readonly(
        "operators",
        &fs::py_flow_graph::operators,
        "Returns the graph's sequence (list) of :class:`FlowOperator` instances.");
    pyfgraph.def_property_readonly("single_flow",
                                   &fs::py_flow_graph::single_flow,
                                   "Returns True if the graph has single direction flow (tree).");

    pyfgraph.def("impl",
                 &fs::py_flow_graph::impl,
                 py::return_value_policy::reference_internal,
                 R"doc(impl(self) -> FlowGraphImpl

        Graph implementation getter.

        Returns
        -------
        impl : :class:`FlowGraphImpl`
            The flow graph implementation instance.

        )doc");

    pyfgraph.def_property_readonly("graph_snapshot_keys",
                                   &fs::py_flow_graph::graph_snapshot_keys,
                                   "Returns the list of graph snasphot names");
    pyfgraph.def("graph_snapshot",
                 &fs::py_flow_graph::graph_snapshot,
                 R"doc(graph_snapshot(self, name: str) -> FlowGraph

        Graph snapshot getter.

        )doc");
    pyfgraph.def_property_readonly("elevation_snapshot_keys",
                                   &fs::py_flow_graph::elevation_snapshot_keys,
                                   "Returns the list of elevation snapshot names");
    pyfgraph.def("elevation_snapshot",
                 &fs::py_flow_graph::elevation_snapshot,
                 R"doc(elevation_snapshot(self, name: str) -> numpy.ndarray

        Elevation snapshot getter.

        The returned array has the same shape than the arrays of the grid used
        to build the flow graph.

        )doc");

    pyfgraph.def("update_routes",
                 &fs::py_flow_graph::update_routes,
                 R"doc(update_routes(self, elevation: numpy.ndarray) -> numpy.ndarray

        (Re-)Compute flow routing and update the flow graph from an input
        topgraphic surface.

        This applies in chain the flow operators and takes snapshots (if any).

        Parameters
        ----------
        elevation : numpy.ndarray
            The input topographic surface used to perform flow routing. The shape of
            the array must match the shape of the arrays of the grid used to build
            the flow graph.

        Returns
        -------
        elevation_copy : numpy.ndarray
            A copy of the input elevation, maybe updated (e.g., with filled depressions).

        )doc");

    pyfgraph.def_property("base_levels",
                          &fs::py_flow_graph::base_levels,
                          &fs::py_flow_graph::set_base_levels,
                          "Indices of the base level nodes.");

    pyfgraph.def_property("mask",
                          &fs::py_flow_graph::mask,
                          &fs::py_flow_graph::set_mask,
                          "Mask where elements with value ``True`` correspond to "
                          "grid nodes that are not included in the flow graph.");

    using data_array_type = fs::py_flow_graph::data_array_type;
    using data_type = fs::py_flow_graph::data_type;

    pyfgraph
        .def("accumulate",
             py::overload_cast<data_array_type&, const data_array_type&>(
                 &fs::py_flow_graph::accumulate, py::const_),
             py::arg("acc"),
             py::arg("src"),
             R"doc(accumulate(*args, **kwargs) -> numpy.ndarray | None

             Traverse the flow graph in the top->down direction and accumulate
             locally produced quantities or fluxes.

             The local quantitites or fluxes (i.e., source) may for example
             correspond to precipitation, surface water runoff, sediment flux,
             etc.

             The accumulated values represent at each node of the graph the
             integral of the source over the node upslope contributing area.

             For example, if the source units are meters (height), the units of
             the output accumulated values are cubic meters (volume).

             Note: the shape of the input and/or output arrays must match the
             shape of the node arrays of the grid used to build the flow graph.

             Overloaded method that supports the following signatures:

             1. ``accumulate(acc: numpy.ndarray, src: numpy.ndarray) -> None``

             2. ``accumulate(acc: numpy.ndarray, src: float) -> None``

             3. ``accumulate(src: numpy.ndarray) -> numpy.ndarray``

             4. ``accumulate(src: float) -> numpy.ndarray``

             Parameters
             ----------
             acc : numpy.ndarray
                 Output accumulated values (array values are overwritten).
                 Only for overloads 1/2.
             src : numpy.ndarray or float
                 Spatially variable (overloads 1/3) or uniform (overloads 2/4)
                 source term (values must be given per area unit).

             Returns
             -------
             acc : numpy.ndarray
                 Output accumulated values, only for overloads 3/4. Otherwise return None.

             )doc")
        .def("accumulate",
             py::overload_cast<data_array_type&, data_type>(&fs::py_flow_graph::accumulate,
                                                            py::const_))
        .def("accumulate",
             py::overload_cast<const data_array_type&>(&fs::py_flow_graph::accumulate, py::const_))
        .def("accumulate",
             py::overload_cast<data_type>(&fs::py_flow_graph::accumulate, py::const_));

    pyfgraph.def("basins",
                 &fs::py_flow_graph::basins,
                 R"doc(basins(self) -> numpy.ndarray

        Delineate catchments (or basins).

        A basin is defined by all adjacent nodes that flow towards the
        same outlet (or pit) graph node.

        Returns
        -------
        basins : numpy.ndarray
            Basin ids. The shape of the array is the same than the shape
            of the arrays of the grid used to build the flow graph.

        Notes
        -----
        Results may be cached so the same computation is not run multiple times.

        All masked grid nodes have the same assigned catchment id set by the maximum
        limit of the integer value range. It is therefore preferable to mask the
        results prior to, e.g., plotting it.

        )doc");

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
