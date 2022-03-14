#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "grid.hpp"
#include "sink_resolver.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


void
add_sink_resolvers_bindings(py::module& m)
{
    fs::detail::register_sink_resolvers<fs::py_profile_grid>();
    fs::detail::register_sink_resolvers<fs::py_raster_grid>();

    // ==== Binding of the SinkResolverMethods enumeration ==== /
    py::enum_<fs::sink_resolver_methods> methods(
        m, "SinkResolverMethods", py::arithmetic(), "Type of flow router algorithm.");
    methods.value("NONE", fs::sink_resolver_methods::none)
        .value("FILL_PFLOOD", fs::sink_resolver_methods::fill_pflood)
        .value("FILL_MST_KRUSKAL", fs::sink_resolver_methods::fill_mst_kruskal)
        .value("FILL_MST_BORUVKA", fs::sink_resolver_methods::fill_mst_boruvka)
        .value("FILL_AUTO", fs::sink_resolver_methods::fill_auto)
        .value("CARVE_MST_KRUSKAL", fs::sink_resolver_methods::carve_mst_kruskal)
        .value("CARVE_MST_BORUVKA", fs::sink_resolver_methods::carve_mst_boruvka);

    // ==== Binding of the BaseSinkResolver class ==== //
    py::class_<fs::detail::sink_resolver_method, std::shared_ptr<fs::detail::sink_resolver_method>>(
        m, "BaseSinkResolver");

    // ==== Binding of the NoSinkResolver class ==== //
    py::class_<fs::detail::no_sink_resolver_method,
               fs::detail::sink_resolver_method,
               std::shared_ptr<fs::detail::no_sink_resolver_method>>(m, "NoSinkResolver")
        .def(py::init());
}
