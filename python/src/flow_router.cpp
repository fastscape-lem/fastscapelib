#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "fastscapelib/profile_grid.hpp"
#include "fastscapelib/raster_grid.hpp"
#include "fastscapelib/flow_router.hpp"

#include "flow_router.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


void
add_flow_routers_bindings(py::module& m)
{
    using namespace fs::detail;

    fs::detail::register_flow_routers<fs::profile_grid>();
    fs::detail::register_flow_routers<fs::raster_grid>();

    // ==== Binding of the FlowRouterMethods enumeration ==== /
    py::enum_<fs::flow_router_methods> methods(
        m, "FlowRouterMethods", py::arithmetic(), "Type of flow router algorithm.");
    methods.value("ONE_CHANNEL", fs::flow_router_methods::one_channel)
        .value("SINGLE", fs::flow_router_methods::single)
        .value("MULTIPLE", fs::flow_router_methods::multiple)
        .value("SINGLE_PARALLEL", fs::flow_router_methods::single_parallel)
        .value("DUMMY", fs::flow_router_methods::dummy);

    // ==== Binding of the BaseFlowRouter class ==== //
    py::class_<flow_router_method, std::shared_ptr<flow_router_method>>(m, "BaseFlowRouter");

    // ==== Binding of the DummyFlowRouter class ==== //
    py::class_<dummy_flow_router_method,
               flow_router_method,
               std::shared_ptr<dummy_flow_router_method>>(m, "DummyFlowRouter")
        .def(py::init());

    // ==== Binding of the SingleFlowRouter class ==== //
    py::class_<single_flow_router_method,
               flow_router_method,
               std::shared_ptr<single_flow_router_method>>(m, "SingleFlowRouter")
        .def(py::init());

    // ==== Binding of the MultipleFlowRouter class ==== //
    py::class_<multiple_flow_router_method,
               flow_router_method,
               std::shared_ptr<multiple_flow_router_method>>(m, "MultipleFlowRouter")
        .def(py::init<double, double>(), py::arg("param1"), py::arg("param2"));
}
