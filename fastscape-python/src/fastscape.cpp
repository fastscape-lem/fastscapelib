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
void compute_receivers_d8_py(xt::pytensor<size_t, 1>& receivers,
                             xt::pytensor<size_t, 1>& dist2receivers,
                             xt::pytensor<T, 2>& elevation,
                             xt::pytensor<bool, 2>& active_nodes,
                             double dx,
                             double dy) {
    py::gil_scoped_release release;
    fs::compute_receivers_d8(receivers, dist2receivers,
                             elevation, active_nodes,
                             dx, dy);
}

template<class T>
void test_py(xt::pytensor<T, 1>& receivers) {
    receivers(438) = 1;
}


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

    m.def("compute_receivers_d8_d", &compute_receivers_d8_py<double>,
          "Compute D8 flow receivers, a single receiver for each grid node.");
    m.def("fill_sinks_flat_d", &fill_sinks_flat_py<double>,
          "Fill depressions in elevation data (flat surfaces).");
    m.def("fill_sinks_sloped_d", &fill_sinks_sloped_py<double>,
          "Fill depressions in elevation data (slightly sloped surfaces).");

    m.def("test_d", &test_py<double>,
          "Compute D8 flow receivers, a single receiver for each grid node.");
}
