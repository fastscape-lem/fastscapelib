#include <memory>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "fastscapelib/flow/sink_resolver.hpp"

#include "grid.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


void
add_sink_resolvers_bindings(py::module& m)
{
    py::class_<fs::no_sink_resolver>(m, "NoSinkResolver").def(py::init());

    py::class_<fs::pflood_sink_resolver>(m, "PFloodSinkResolver").def(py::init());


    py::enum_<fs::mst_method> mst_method(m, "MSTMethod", py::arithmetic());
    mst_method.value("KRUSKAL", fs::mst_method::kruskal).value("BORUVKA", fs::mst_method::boruvka);

    py::enum_<fs::mst_route_method> mst_route_method(m, "MSTRouteMethod", py::arithmetic());
    mst_route_method.value("BASIC", fs::mst_route_method::basic)
        .value("CARVE", fs::mst_route_method::carve);

    py::class_<fs::mst_sink_resolver>(m, "MSTSinkResolver")
        .def(py::init())
        .def(py::init<fs::mst_method, fs::mst_route_method>())
        .def_readwrite("basin_method", &fs::mst_sink_resolver::basin_method)
        .def_readwrite("route_method", &fs::mst_sink_resolver::route_method);
}
