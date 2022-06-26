#include <memory>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "fastscapelib/flow/flow_router.hpp"

#include "grid.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


void
add_flow_routers_bindings(py::module& m)
{
    py::class_<fs::single_flow_router>(m, "SingleFlowRouter").def(py::init());

    py::class_<fs::multiple_flow_router>(m, "MultipleFlowRouter")
        .def(py::init<double, double>())
        .def_readwrite("p1", &fs::multiple_flow_router::p1)
        .def_readwrite("p2", &fs::multiple_flow_router::p2);
}
