/**
 * @file
 * @brief Fastscapelib Python bindings.
 */
#include "pybind11/pybind11.h"
#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"

#include "fastscapelib/version.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


void
add_grid_bindings(py::module_&);
void
add_flow_graph_bindings(py::module_&);
void
add_eroders_bindings(py::module_&);


PYBIND11_MODULE(_fastscapelib_py, m)
{
    m.doc() = "A collection of efficient algorithms"
              "for processing topographic data and landscape evolution modeling.";

    xt::import_numpy();

    m.attr("__version__") = fs::version::version_str;

    py::module grid_m = m.def_submodule("grid", "Fastscapelib's grid module");
    add_grid_bindings(grid_m);

    py::module flow_m = m.def_submodule("flow", "Fastscapelib's flow routing module");
    add_flow_graph_bindings(flow_m);

    py::module eroders_m = m.def_submodule("eroders", "Fastscapelib's erosion module");
    add_eroders_bindings(eroders_m);
}
