#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "xtensor-python/pytensor.hpp"

#include "fastscapelib/algo/pflood.hpp"

#include "grid.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


void
add_sinks_bindings(py::module& m)
{
    using pytensor_1d = xt::pytensor<double, 1>;
    using pytensor_2d = xt::pytensor<double, 2>;

    m.def("fill_sinks_flat",
          &fs::detail::fill_sinks_flat_impl<fs::py_raster_grid, pytensor_2d>,
          py::call_guard<py::gil_scoped_release>(),
          "Fill depressions in 2D elevation data (flat surfaces).");

    m.def("fill_sinks_flat",
          &fs::detail::fill_sinks_flat_impl<fs::py_profile_grid, pytensor_1d>,
          py::call_guard<py::gil_scoped_release>(),
          "Fill depressions in 1D elevation data (flat surfaces).");

    m.def("fill_sinks_sloped",
          &fs::detail::fill_sinks_sloped_impl<fs::py_raster_grid, pytensor_2d>,
          py::call_guard<py::gil_scoped_release>(),
          "Fill depressions in 2D elevation data (slightly sloped surfaces).");

    m.def("fill_sinks_sloped",
          &fs::detail::fill_sinks_sloped_impl<fs::py_profile_grid, pytensor_1d>,
          py::call_guard<py::gil_scoped_release>(),
          "Fill depressions in 1D elevation data (slightly sloped surfaces).");
}
