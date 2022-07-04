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
}
