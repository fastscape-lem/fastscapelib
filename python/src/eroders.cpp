#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "fastscapelib/eroders/spl.hpp"
#include "fastscapelib/eroders/diffusion_adi.hpp"
#include "fastscapelib/grid/raster_grid.hpp"

#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"

#include "flow_graph.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


template <class K, class T>
auto
erode_stream_power_py(xt::pyarray<T, xt::layout_type::row_major>& erosion,
                      const xt::pyarray<T, xt::layout_type::row_major>& elevation,
                      const xt::pyarray<T, xt::layout_type::row_major>& drainage_area,
                      fs::py_flow_graph& flow_graph,
                      const K k_coef,
                      double m_exp,
                      double n_exp,
                      double dt,
                      double tolerance) -> std::size_t
{
    py::gil_scoped_release release;
    return fs::erode_stream_power(
        erosion, elevation, drainage_area, flow_graph, k_coef, m_exp, n_exp, dt, tolerance);
}

template <class K, class T>
void
erode_linear_diffusion_py(xt::pytensor<T, 2>& erosion,
                          const xt::pytensor<T, 2>& elevation,
                          const K k_coef,
                          double dt,
                          double dx,
                          double dy)
{
    py::gil_scoped_release release;
    fs::erode_linear_diffusion(erosion, elevation, k_coef, dt, dx, dy);
}


void
add_spl_bindings(py::module& m)
{
    m.def("erode_stream_power_d",
          &erode_stream_power_py<double, double>,
          "Compute bedrock channel erosion during a single time step "
          "using the Stream Power Law.");

    m.def("erode_stream_power_var_d",
          &erode_stream_power_py<xt::pyarray<double, xt::layout_type::row_major>, double>,
          "Compute bedrock channel erosion during a single time step "
          "using the Stream Power Law.\n\n"
          "Version with spatially variable stream power coefficient.");
}


void
add_diffusion_adi_bindings(py::module& m)
{
    m.def("erode_linear_diffusion_d",
          &erode_linear_diffusion_py<double, double>,
          "Compute hillslope erosion by linear diffusion on a 2-d regular "
          "grid using finite differences with an Alternating Direction"
          "Implicit (ADI) scheme.");

    m.def("erode_linear_diffusion_var_d",
          &erode_linear_diffusion_py<xt::pytensor<double, 2>&, double>,
          "Compute hillslope erosion by linear diffusion on a 2-d regular "
          "grid using finite differences with an Alternating Direction"
          "Implicit (ADI) scheme.\n\n"
          "Version with spatially variable diffusion coefficient.");
}


void
add_eroders_bindings(py::module& m)
{
    add_spl_bindings(m);
    add_diffusion_adi_bindings(m);
}
