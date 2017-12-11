/**
 * @file
 * @brief Fastscape Python bindings.
*/
#include <cstdint>

#include "pybind11/pybind11.h"
#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pytensor.hpp"

#include "fastscape/utils.hpp"
#include "fastscape/fastscape.hpp"


namespace py = pybind11;
namespace fs = fastscape;


template<class T>
void compute_receivers_d8_py(xt::pytensor<index_t, 1>& receivers,
                             xt::pytensor<T, 1>& dist2receivers,
                             const xt::pytensor<T, 2>& elevation,
                             const xt::pytensor<bool, 2>& active_nodes,
                             double dx,
                             double dy) {
    py::gil_scoped_release release;
    fs::compute_receivers_d8(receivers, dist2receivers,
                             elevation, active_nodes,
                             dx, dy);
}


void compute_donors_py(xt::pytensor<index_t, 1>& ndonors,
                       xt::pytensor<index_t, 2>& donors,
                       const xt::pytensor<index_t, 1>& receivers) {
    py::gil_scoped_release release;
    fs::compute_donors(ndonors, donors, receivers);
}


void compute_stack_py(xt::pytensor<index_t, 1>& stack,
                      const xt::pytensor<index_t, 1>& ndonors,
                      const xt::pytensor<index_t, 2>& donors,
                      const xt::pytensor<index_t, 1>& receivers) {
    py::gil_scoped_release release;
    fs::compute_stack(stack, ndonors, donors, receivers);
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
        "for processing topographic data and landscape evolution modeling.";

    xt::import_numpy();

    m.def("compute_receivers_d8_d", &compute_receivers_d8_py<double>,
          "Compute D8 flow receivers, a single receiver for each grid node.");

    m.def("compute_donors", &compute_donors_py,
          "Compute flow donors (invert flow receivers).");

    m.def("compute_stack", &compute_stack_py,
          "Build a stack of grid node ids such that each node in the stack "
          "is visited before any of its upstream nodes.");

    m.def("fill_sinks_flat_d", &fill_sinks_flat_py<double>,
          "Fill depressions in elevation data (flat surfaces).");

    m.def("fill_sinks_sloped_d", &fill_sinks_sloped_py<double>,
          "Fill depressions in elevation data (slightly sloped surfaces).");
}
