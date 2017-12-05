/*
  <%
  cfg['dependencies'] = ['include/fastscape/sinks.hpp', 'include/fastscape/consts.hpp']
  cfg['compiler_args'] = ['-std=c++14']
  cfg['include_dirs'] = ['/Users/bbovy/miniconda3/envs/fastscape_py36/lib/python3.6/site-packages/numpy/core/include']
  setup_pybind11(cfg)
  %>
*/
#include "pybind11/pybind11.h"
#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pytensor.hpp"

#include "include/fastscape/sinks.hpp"


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


PYBIND11_MODULE(pyfastscape, m) {
    m.doc() = "A collection of efficient algorithms"
        "for processing topographic data and landscape evolution modeling";

    xt::import_numpy();

    m.def("fill_sinks_flat_d", &fill_sinks_flat_py<double>,
          "Fill depressions in elevation data (flat surfaces).");
    m.def("fill_sinks_sloped_d", &fill_sinks_sloped_py<double>,
          "Fill depressions in elevation data (slightly sloped surfaces).");
}
