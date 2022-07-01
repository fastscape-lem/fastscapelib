#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "fastscapelib/eroders/spl.hpp"
#include "fastscapelib/eroders/diffusion_adi.hpp"
#include "fastscapelib/grid/raster_grid.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"

#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"

#include "grid.hpp"
#include "flow_graph.hpp"
#include "pytensor_utils.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


void
add_spl_bindings(py::module& m)
{
    using py_spl_eroder = fs::spl_eroder<fs::py_flow_graph, fs::py_selector>;
    using data_array_type = py_spl_eroder::data_array_type;

    py::class_<py_spl_eroder>(m, "SPLEroder")
        .def(py::init<fs::py_flow_graph&, double, double, double, double>())
        .def(py::init<fs::py_flow_graph&, data_array_type&, double, double, double>())
        .def_property("k_coef",
                      &py_spl_eroder::k_coef,
                      [](py_spl_eroder& self, py::object value)
                      {
                          if (py::isinstance<py::float_>(value))
                          {
                              self.set_k_coef(value.cast<double>());
                          }
                          else if (py::isinstance<data_array_type>(value))
                          {
                              self.set_k_coef(value.cast<data_array_type>());
                          }
                      })
        .def_property("area_exp", &py_spl_eroder::area_exp, &py_spl_eroder::set_area_exp)
        .def_property("slope_exp", &py_spl_eroder::slope_exp, &py_spl_eroder::set_slope_exp)
        .def_property_readonly("tolerance", &py_spl_eroder::tolerance)
        .def_property_readonly("n_corr", &py_spl_eroder::n_corr)
        .def("erode", &py_spl_eroder::erode, py::call_guard<py::gil_scoped_release>());
}


void
add_diffusion_adi_bindings(py::module& m)
{
    using py_diffusion_adi_eroder = fs::diffusion_adi_eroder<fs::py_raster_grid, fs::py_selector>;
    using data_array_type = py_diffusion_adi_eroder::data_array_type;

    py::class_<py_diffusion_adi_eroder>(m, "DiffusionADIEroder")
        .def(py::init<fs::py_raster_grid&, double>())
        .def(py::init<fs::py_raster_grid&, data_array_type&>())
        .def_property("k_coef",
                      &py_diffusion_adi_eroder::k_coef,
                      [](py_diffusion_adi_eroder& self, py::object value)
                      {
                          if (py::isinstance<py::float_>(value))
                          {
                              self.set_k_coef(value.cast<double>());
                          }
                          else if (py::isinstance<data_array_type>(value))
                          {
                              self.set_k_coef(value.cast<data_array_type>());
                          }
                      })
        .def("erode", &py_diffusion_adi_eroder::erode, py::call_guard<py::gil_scoped_release>());
}


void
add_eroders_bindings(py::module& m)
{
    add_spl_bindings(m);
    add_diffusion_adi_bindings(m);
}
