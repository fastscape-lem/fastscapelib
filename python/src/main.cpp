/**
 * @file
 * @brief Fastscapelib Python bindings.
*/
#include <cstdint>

#include "pybind11/pybind11.h"
#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"

#include "fastscapelib/utils.hpp"
#include "fastscapelib/fastscapelib.hpp"

#include "grid.cpp"
#include "modelings.cpp"


namespace py = pybind11;
namespace fs = fastscapelib;


PYBIND11_MODULE(_fastscapelib_py, m)
{
    m.doc() = "A collection of efficient algorithms"
        "for processing topographic data and landscape evolution modeling.";

    xt::import_numpy();

    m.attr("__version__") = fs::version::version_str;

    add_grid_bindings(m);
    add_modelings_bindings(m);
}
