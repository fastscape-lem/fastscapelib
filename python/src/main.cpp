/**
 * @file
 * @brief Fastscapelib Python bindings.
 */
#include <cstdint>

#include "pybind11/pybind11.h"
#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"

#include "fastscapelib/utils.hpp"
#include "fastscapelib/fastscapelib.hpp"

#include "grid.cpp"
#include "modelings.cpp"
#include "flow_router.cpp"
#include "sink_resolver.cpp"
#include "sinks.cpp"
#include "flow_graph.cpp"


namespace py = pybind11;
namespace fs = fastscapelib;


PYBIND11_MODULE(_fastscapelib_py, m)
{
    m.doc() = "A collection of efficient algorithms"
              "for processing topographic data and landscape evolution modeling.";

    xt::import_numpy();

    m.attr("__version__") = fs::version::version_str;

    py::module grid_m = m.def_submodule("grid", "The grid module of Fastscapelib");
    add_grid_bindings(grid_m);

    py::module flow_graph_m
        = m.def_submodule("flow_graph", "The flow graph module of Fastscapelib");
    add_flow_routers_bindings(flow_graph_m);
    add_sink_resolvers_bindings(flow_graph_m);
    add_flow_graph_bindings(flow_graph_m);

    // TODO: merge flow_graph and sinks modules?
    py::module sinks_m = m.def_submodule("sinks", "Various algorithms for sink filling");
    add_sinks_bindings(sinks_m);

    py::module algo_m = m.def_submodule("algo", "The algorithm module of Fastscapelib");
    add_modelings_bindings(algo_m);
}
