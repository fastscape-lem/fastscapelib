#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "fastscapelib/flow/flow_router.hpp"

#include "grid.hpp"
#include "flow_router.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


void
add_flow_routers_bindings(py::module& m)
{
    py::class_<fs::detail::py_flow_router>(m, "FlowRouter");

    py::class_<fs::detail::py_dummy_flow_router, fs::detail::py_flow_router>(m, "DummyFlowRouter")
        .def(py::init());

    py::class_<fs::detail::py_single_flow_router, fs::detail::py_flow_router>(m, "SingleFlowRouter")
        .def(py::init());

    py::class_<fs::detail::py_multiple_flow_router, fs::detail::py_flow_router>(
        m, "MultipleFlowRouter")
        .def(py::init<double, double>(), py::arg("param1"), py::arg("param2"));
}
