/*
  <%
  cfg['dependencies'] = ['depressions.hpp', 'constants.hpp']
  cfg['compiler_args'] = ['-std=c++14']
  cfg['include_dirs'] = ['/Users/bbovy/miniconda3/envs/fastscape_py36/lib/python3.6/site-packages/numpy/core/include']
  setup_pybind11(cfg)
  %>
*/
#include "pybind11/pybind11.h"
#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pyarray.hpp"

#include "depressions.hpp"


namespace py = pybind11;


template<class T>
void priority_flood_original_py(xt::pyarray<T>& elevation) {
    py::gil_scoped_release release;
    priority_flood_original(elevation);
}


PYBIND11_MODULE(pyfastscape, m) {
    m.doc() = "A collection of efficient algorithms for"
        "processing topographic data and landscape evolution modeling";

    xt::import_numpy();

    m.def("priority_flood_original_d", &priority_flood_original_py<double>,
          "Fill depressions in elevation data (flat surfaces).");
    //m.def("priority_flood_epsilon_d", &priority_flood_epsilon_py<double>,
    //      "Fill depressions in elevation data (no flat surfaces)");
}
