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
    py::options options;
    options.disable_function_signatures();

    using py_spl_eroder = fs::spl_eroder<fs::py_flow_graph, fs::py_selector>;
    using data_array_type = py_spl_eroder::data_array_type;

    py::class_<py_spl_eroder> spl_eroder(
        m,
        "SPLEroder",
        R"doc(Bedrock channel erosion modelled using the Stream Power Law.

        It numerically solves the Stream Power Law (SPL) using an implicit finite
        difference scheme 1st order in space and time. The method is detailed in
        Braun and Willet's (2013) and has been slightly adapted.

        For more details about the Stream Power Law and the numerical scheme used
        here, see :cpp:class:`~fastscapelib::spl_eroder` from the C++ API. See also
        :cite:t:`Braun2013`.

        )doc");

    spl_eroder.def(
        py::init<fs::py_flow_graph&, double, double, double, double>(),
        py::arg("flow_graph"),
        py::arg("k_coef"),
        py::arg("area_exp"),
        py::arg("slope_exp"),
        py::arg("tolerance") = 1e-3,
        R"doc(__init__(self, flow_graph: FlowGraph, k_coef: Union[float, numpy.ndarray], area_exp: float, slope_exp: float, tolerance: float = 1e-3)

            SPLEroder initializer.

            Parameters
            ----------
            flow_graph : :class:`FlowGraph`
                Flow graph instance.
            k_coef : float or numpy.ndarray
                Spatially uniform or variable erosion coefficient value. For the
                spatially variable case, the shape of the array must match the
                shape of the grid used (node arrays) or the size of the graph.
            area_exp : float
                Drainage area exponent value.
            slope_exp : float
                Channel slope exponent value.
            tolerance : float, optional
                Tolerance controlling the convergence of the Newton-Raphson
                iterations used for the non-linear case (i.e., ``slope_exp != 1``).

            )doc");
    spl_eroder.def(py::init<fs::py_flow_graph&, data_array_type&, double, double, double>(),
                   py::arg("flow_graph"),
                   py::arg("k_coef"),
                   py::arg("area_exp"),
                   py::arg("slope_exp"),
                   py::arg("tolerance") = 1e-3);

    spl_eroder
        .def_property(
            "k_coef",
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
            },
            "Erosion coefficient (spatially variable or uniform).")
        .def_property("area_exp",
                      &py_spl_eroder::area_exp,
                      &py_spl_eroder::set_area_exp,
                      "Drainage area exponent value.")
        .def_property("slope_exp",
                      &py_spl_eroder::slope_exp,
                      &py_spl_eroder::set_slope_exp,
                      "Channel slope exponent value.")
        .def_property_readonly(
            "tolerance",
            &py_spl_eroder::tolerance,
            R"doc(Tolerance controlling the convergence of the Newton-Raphson method
                               used for the non-linear case (i.e., ``slope_exp != 1)``.)doc")
        .def_property_readonly(
            "n_corr",
            &py_spl_eroder::n_corr,
            R"doc(Returns the number of nodes for which erosion has been arbitrarily
                               limited during the last computed time step
                               (see :cpp:func:`~fastscapelib::spl_eroder::n_corr`) from the C++ API.)doc");

    spl_eroder.def(
        "erode",
        &py_spl_eroder::erode,
        py::call_guard<py::gil_scoped_release>(),
        py::arg("elevation"),
        py::arg("drainage_area"),
        py::arg("dt"),
        R"doc(erode(elevation: numpy.ndarray, drainage_area: numpy.ndarray, dt: float) -> numpy.ndarray

            Slove SPL for one time step.

            Parameters
            ----------
            elevation : numpy.ndarray
                Surface topography elevation at each grid node (shape must match
                grid node array shape).
            drainage_area : numpy.ndarray
                Upslope drainage area at each grid node (same shape than ``elevation``).
            dt : float
                Duration of the time step.

            Returns
            -------
            erosion : numpy.ndarray
                SPL vertical erosion computed for the time step (>= 0). Same shape than
                the input arrays.

            )doc");
}


void
add_diffusion_adi_bindings(py::module& m)
{
    py::options options;
    options.disable_function_signatures();

    using py_diffusion_adi_eroder = fs::diffusion_adi_eroder<fs::py_raster_grid, fs::py_selector>;
    using data_array_type = py_diffusion_adi_eroder::data_array_type;

    py::class_<py_diffusion_adi_eroder> diffusion_adi_eroder(
        m,
        "DiffusionADIEroder",
        R"doc(Hillslope erosion using linear diffusion.

        It numerically solves the diffusion equation using an Alternating
        Direction Implicit (ADI) scheme.

        The equation is given by:

        .. math::

           \frac{\partial h}{\partial t} = K \nabla^2 h

        where :math:`K` is an erosion coefficient (diffusivity) and :math:`\nabla^2 h`
        is the local curvature of the topographic surface.

        This equation implies that the amount of sediment eroded is linearly
        proportional to the local gradient of the topographic surface.

        Notes: only raster grids are supported. This eroder assumes Dirichlet
        boundary conditions at the border nodes of the raster grid.

        )doc");

    diffusion_adi_eroder.def(
        py::init<fs::py_raster_grid&, double>(),
        py::arg("grid"),
        py::arg("k_coef"),
        R"doc(__init__(self, grid: RasterGrid, k_coef: Union[float, numpy.ndarray]) -> None

        DiffusionADIEroder initializer.

        Parameters
        ----------
        grid : :class:`RasterGrid`
            A raster grid instance.
        k_coef : float or numpy.ndarray
            Spatially uniform or variable erosion coefficient (diffusivity) value.
            For the spatially variable case, the shape of the array must match the
            shape of the raster grid.

        )doc");
    diffusion_adi_eroder.def(
        py::init<fs::py_raster_grid&, data_array_type&>(), py::arg("grid"), py::arg("k_coef"));

    diffusion_adi_eroder.def_property(
        "k_coef",
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
        },
        "Erosion coefficient (spatially variable or uniform).");

    diffusion_adi_eroder.def("erode",
                             &py_diffusion_adi_eroder::erode,
                             py::call_guard<py::gil_scoped_release>(),
                             py::arg("elevation"),
                             py::arg("dt"),
                             R"doc(erode(elevation: numpy.ndarray, dt: float) -> numpy.ndarray

            Slove diffusion for one time step.

            Parameters
            ----------
            elevation : numpy.ndarray
                Surface topography elevation at each grid node (shape must match
                the shape of the raster grid).
            dt : float
                Duration of the time step.

            Returns
            -------
            erosion : numpy.ndarray
                Diffusion vertical erosion (> 0) or accumulation (< 0) computed
                for the time step. Same shape than the input ``elevation`` array.

            )doc");
}


void
add_eroders_bindings(py::module& m)
{
    add_spl_bindings(m);
    add_diffusion_adi_bindings(m);
}
