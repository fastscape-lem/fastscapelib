/**
 * @file
 * @brief Fastscapelib Python bindings.
*/
#include <cstdint>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"

#include "fastscapelib/utils.hpp"
#include "fastscapelib/fastscapelib.hpp"

#include "xtensor_utils.hpp"


namespace py = pybind11;
using namespace pybind11::literals;  // use the `_a` literal

namespace fs = fastscapelib;


py::dict get_versions()
{
    return py::dict("version"_a=fs::version::version,
                    "git_hash_full"_a=fs::version::git_hash_full);
}


template<class T>
void fill_sinks_flat_py(xt::pytensor<T, 2>& elevation)
{
    py::gil_scoped_release release;
    fs::fill_sinks_flat(elevation);
}


template<class T>
void fill_sinks_sloped_py(xt::pytensor<T, 2>& elevation)
{
    py::gil_scoped_release release;
    fs::fill_sinks_sloped(elevation);
}


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


PYBIND11_MODULE(_fastscapelib_py, m)
{
    m.doc() = "A collection of efficient algorithms"
        "for processing topographic data and landscape evolution modeling.";

    xt::import_numpy();

    m.def("get_versions", &get_versions,
          "Get version info.");

    //*******
    //* Grid
    //*******

    using profile_grid_py = fs::profile_grid_xt<fs::pytensor_selector>;

    py::enum_<fs::node_status>(m, "NodeStatus", py::arithmetic(),
                               "Status of grid/mesh nodes either inside the domain "
                               "or on the domain boundary.")
        .value("CORE", fs::node_status::core)
        .value("FIXED_VALUE_BOUNDARY", fs::node_status::fixed_value_boundary)
        .value("FIXED_GRADIENT_BOUNDARY", fs::node_status::fixed_gradient_boundary)
        .value("LOOPED_BOUNDARY", fs::node_status::looped_boundary);

    py::class_<fs::edge_nodes_status>(m, "EdgeNodesStatus")
        .def(py::init<fs::node_status, fs::node_status>())
        .def_readonly("left", &fs::edge_nodes_status::left)
        .def_readonly("right", &fs::edge_nodes_status::right);

    py::class_<fs::node>(m, "Node")
        .def(py::init<std::size_t, fs::node_status>())
        .def_readonly("idx", &fs::node::idx)
        .def_readonly("status", &fs::node::status);

    py::class_<profile_grid_py>(m, "ProfileGrid")
        .def(py::init<std::size_t, double, const fs::edge_nodes_status,
                      const std::vector<fs::node>&>())
        .def_property_readonly("size", &profile_grid_py::size)
        .def_property_readonly("spacing", &profile_grid_py::spacing)
        .def_property_readonly("node_status", &profile_grid_py::node_status);

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

    m.def("fill_sinks_flat_d", &fill_sinks_flat_py<double>,
          "Fill depressions in elevation data (flat surfaces).");

    m.def("fill_sinks_sloped_d", &fill_sinks_sloped_py<double>,
          "Fill depressions in elevation data (slightly sloped surfaces).");

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
