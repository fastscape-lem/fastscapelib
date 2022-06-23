#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "grid.hpp"
#include "sink_resolver.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


void
add_sink_resolvers_bindings(py::module& m)
{
    py::class_<fs::detail::py_sink_resolver>(m, "SinkResolver");

    py::class_<fs::detail::py_no_sink_resolver, fs::detail::py_sink_resolver>(m, "NoSinkResolver")
        .def(py::init());
}
