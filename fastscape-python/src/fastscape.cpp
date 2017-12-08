/**
 * @file
 * @brief Fastscape Python bindings.
*/
#include "pybind11/pybind11.h"
#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pytensor.hpp"

#include "fastscape/fastscape.hpp"


namespace py = pybind11;
namespace fs = fastscape;


template<class T>
void fill_sinks_flat_py(xt::pytensor<T, 2>& elevation) {
    py::gil_scoped_release release;
    fs::fill_sinks_flat(elevation);
}


template<class T>
void fill_sinks_sloped_py(xt::pytensor<T, 2>& elevation) {
    py::gil_scoped_release release;
    fs::fill_sinks_sloped(elevation);
}


PYBIND11_MODULE(fastscape, m) {
    m.doc() = "A collection of efficient algorithms"
        "for processing topographic data and landscape evolution modeling";

    xt::import_numpy();

    m.def("fill_sinks_flat_d", &fill_sinks_flat_py<double>,
          "Fill depressions in elevation data (flat surfaces).");
    m.def("fill_sinks_sloped_d", &fill_sinks_sloped_py<double>,
          "Fill depressions in elevation data (slightly sloped surfaces).");
}
