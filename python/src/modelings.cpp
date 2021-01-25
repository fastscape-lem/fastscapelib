/**
 * @file
 * @brief Fastscapelib modelings Python bindings.
 */
#include "fastscapelib/bedrock_channel.hpp"
#include "fastscapelib/hillslope.hpp"
#include "fastscapelib/flow_routing.hpp"

#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

namespace py = pybind11;
namespace fs = fastscapelib;


template<class T>
void compute_receivers_d8_py(xt::pytensor<index_t, 1>& receivers,
                             xt::pytensor<T, 1>& dist2receivers,
                             const xt::pytensor<T, 2>& elevation,
                             const xt::pytensor<bool, 2>& active_nodes,
                             double dx,
                             double dy)
{
    py::gil_scoped_release release;
    fs::compute_receivers_d8(receivers, dist2receivers,
                             elevation, active_nodes,
                             dx, dy);
}


void compute_donors_py(xt::pytensor<index_t, 1>& ndonors,
                       xt::pytensor<index_t, 2>& donors,
                       const xt::pytensor<index_t, 1>& receivers)
{
    py::gil_scoped_release release;
    fs::compute_donors(ndonors, donors, receivers);
}


void compute_stack_py(xt::pytensor<index_t, 1>& stack,
                      const xt::pytensor<index_t, 1>& ndonors,
                      const xt::pytensor<index_t, 2>& donors,
                      const xt::pytensor<index_t, 1>& receivers)
{
    py::gil_scoped_release release;
    fs::compute_stack(stack, ndonors, donors, receivers);
}


index_t compute_basins_py(xt::pytensor<index_t, 1>& basins,
                          xt::pytensor<index_t, 1>& outlets_or_pits,
                          const xt::pytensor<index_t, 1>& stack,
                          const xt::pytensor<index_t, 1>& receivers)
{
    py::gil_scoped_release release;
    return fs::compute_basins(basins, outlets_or_pits, stack, receivers);
}


index_t find_pits_py(xt::pytensor<index_t, 1>& pits,
                     const xt::pytensor<index_t, 1>& outlets_or_pits,
                     const xt::pytensor<bool, 2>& active_nodes,
                     index_t nbasins)
{
    py::gil_scoped_release release;
    return fs::find_pits(pits, outlets_or_pits, active_nodes, nbasins);
}


template<class T>
void compute_drainage_area_mesh_py(xt::pytensor<T, 1>& drainage_area,
                                   const xt::pytensor<T, 1>& cell_area,
                                   const xt::pytensor<index_t, 1>& stack,
                                   const xt::pytensor<index_t, 1>& receivers)
{
    py::gil_scoped_release release;
    fs::compute_drainage_area(drainage_area, cell_area, stack, receivers);
}


template<class T>
void compute_drainage_area_grid_py(xt::pytensor<T, 2>& drainage_area,
                                   const xt::pytensor<index_t, 1>& stack,
                                   const xt::pytensor<index_t, 1>& receivers,
                                   double dx,
                                   double dy)
{
    py::gil_scoped_release release;
    fs::compute_drainage_area(drainage_area, stack, receivers, dx, dy);
}


template<class K, class T>
index_t erode_stream_power_py(xt::pyarray<T>& erosion,
                              const xt::pyarray<T>& elevation,
                              const xt::pytensor<index_t, 1>& stack,
                              const xt::pytensor<index_t, 1>& receivers,
                              const xt::pytensor<T, 1>& dist2receivers,
                              const xt::pyarray<T>& drainage_area,
                              const K k_coef,
                              double m_exp,
                              double n_exp,
                              double dt,
                              double tolerance)
{
    py::gil_scoped_release release;
    return fs::erode_stream_power(erosion, elevation,
                                  stack, receivers,
                                  dist2receivers, drainage_area,
                                  k_coef, m_exp, n_exp,
                                  dt, tolerance);
}


template<class K, class T>
void erode_linear_diffusion_py(xt::pytensor<T, 2>& erosion,
                               const xt::pytensor<T, 2>& elevation,
                               const K k_coef,
                               double dt,
                               double dx,
                               double dy)
{
    py::gil_scoped_release release;
    fs::erode_linear_diffusion(erosion, elevation, k_coef, dt, dx, dy);
}

void add_modelings_bindings(py::module& m)
{
    m.def("compute_receivers_d8_d", &compute_receivers_d8_py<double>,
          "Compute D8 flow receivers, a single receiver for each grid node.");

    m.def("compute_donors", &compute_donors_py,
          "Compute flow donors (invert flow receivers).");

    m.def("compute_stack", &compute_stack_py,
          "Build a stack of grid node ids such that each node in the stack "
          "is visited before any of its upstream nodes.");

    m.def("compute_basins", &compute_basins_py,
          "Compute basin ids. Return total number of basins");

    m.def("find_pits", &find_pits_py,
          "Find pit node ids. Return total number of pit nodes");

    m.def("compute_drainage_area_mesh_d",
          &compute_drainage_area_mesh_py<double>,
          "Compute drainage area on a mesh.");

    m.def("compute_drainage_area_grid_d",
          &compute_drainage_area_grid_py<double>,
          "Compute drainage area on a 2D grid.");

    m.def("erode_stream_power_d", &erode_stream_power_py<double, double>,
          "Compute bedrock channel erosion during a single time step "
          "using the Stream Power Law.");

    m.def("erode_stream_power_var_d", &erode_stream_power_py<xt::pyarray<double>&, double>,
          "Compute bedrock channel erosion during a single time step "
          "using the Stream Power Law.\n\n"
          "Version with spatially variable stream power coefficient.");

    m.def("erode_linear_diffusion_d", &erode_linear_diffusion_py<double, double>,
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
 
