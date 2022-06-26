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
    py::class_<fs::py_flow_graph> pyfgraph(m, "FlowGraph");

    fs::add_init_methods<fs::py_profile_grid>(pyfgraph);
    fs::add_init_methods<fs::py_raster_grid>(pyfgraph);

    pyfgraph.def("update_routes", &fs::py_flow_graph::update_routes)
        .def("receivers", &fs::py_flow_graph::receivers)
        .def("receivers_count", &fs::py_flow_graph::receivers_count)
        .def("receivers_distance", &fs::py_flow_graph::receivers_distance)
        .def("receivers_weight", &fs::py_flow_graph::receivers_weight)
        .def("donors", &fs::py_flow_graph::donors)
        .def("donors_count", &fs::py_flow_graph::donors_count)
        .def("dfs_stack", &fs::py_flow_graph::dfs_stack);

    pyfgraph
        .def("accumulate",
             py::overload_cast<const xt::pyarray<double, xt::layout_type::row_major>&>(
                 &fs::py_flow_graph::accumulate, py::const_))
        .def("accumulate",
             py::overload_cast<const double&>(&fs::py_flow_graph::accumulate, py::const_));
}
