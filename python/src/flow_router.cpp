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

    py::class_<fs::multi_flow_router>(m, "MultiFlowRouter")
        .def(py::init<double>())
        .def_readwrite("slope_exp", &fs::multi_flow_router::slope_exp);

    py::class_<fs::singlemulti_flow_router>(m, "SingleMultiFlowRouter")
        .def(py::init<double>())
        .def_readwrite("slope_exp", &fs::singlemulti_flow_router::slope_exp);
}
