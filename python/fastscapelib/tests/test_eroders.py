import numpy as np
import pytest

from fastscapelib.eroders import DiffusionADIEroder, SPLEroder
from fastscapelib.flow import FlowGraph, MultiFlowRouter, SingleFlowRouter
from fastscapelib.grid import NodeStatus, ProfileGrid, RasterGrid


class TestSPLEroder:
    def test_constructor_properties(self) -> None:
        grid = RasterGrid([2, 2], [1.0, 1.0], NodeStatus.FIXED_VALUE)
        flow_graph = FlowGraph(grid, [SingleFlowRouter()])

        eroder = SPLEroder(flow_graph, 1e-3, 0.4, 1, 1e-5)

        np.testing.assert_equal(eroder.k_coef, np.full(grid.size, 1e-3))
        assert eroder.area_exp == 0.4
        assert eroder.slope_exp == 1
        assert eroder.tolerance == 1e-5

        eroder.area_exp = 0.5
        assert eroder.area_exp == 0.5

        eroder.slope_exp = 1.2
        assert eroder.slope_exp == 1.2

        k_coef_arr = np.arange(grid.size, dtype=np.double).reshape(grid.shape)
        eroder.k_coef = k_coef_arr
        np.testing.assert_equal(eroder.k_coef, k_coef_arr.flatten())

        k_coef = 1e-5
        eroder.k_coef = k_coef
        np.testing.assert_equal(eroder.k_coef, np.full(grid.size, 1e-5))

        with pytest.raises(RuntimeError, match=".*shape mismatch"):
            eroder.k_coef = np.ones((4, 5))

        # multiple flow direction and n != 1 is not supoprted
        flow_graph2 = FlowGraph(grid, [MultiFlowRouter(1.0)])
        with pytest.raises(ValueError, match=".*slope exponent != 1 is not supported"):
            SPLEroder(flow_graph2, 1e-3, 0.4, 1.5, 1e-3)

    @pytest.mark.parametrize("k_coef", [1e-3, np.full((4), 1e-3)])
    def test_profile_grid(self, k_coef) -> None:
        spacing = 300.0
        grid = ProfileGrid(4, 300, NodeStatus.FIXED_VALUE)

        flow_graph = FlowGraph(grid, [SingleFlowRouter()])

        area_exp = 0.5
        slope_exp = 1.0
        tolerance = 1e-3
        eroder = SPLEroder(flow_graph, k_coef, area_exp, slope_exp, tolerance)

        h = 1.0
        elevation = np.array([0.0, h, h, 0.0], dtype="d")

        flow_graph.update_routes(elevation)
        drainage_area = flow_graph.accumulate(1.0)

        dt = 1.0  # use small time step (compare with explicit scheme)
        erosion = eroder.erode(elevation, drainage_area, dt)

        slope = h / spacing
        a = spacing
        k_coef = 1e-3
        err = dt * k_coef * a**area_exp * slope**slope_exp
        expected_erosion = np.array([0.0, err, err, 0.0], dtype="d")

        np.testing.assert_allclose(erosion, expected_erosion, atol=1e-5)
        assert eroder.n_corr == 0

    @pytest.mark.parametrize("k_coef", [1e-3, np.full((2, 2), 1e-3)])
    def test_raster_grid(self, k_coef) -> None:
        # Test on a tiny (2x2) 2-d square grid with a planar surface
        # tilted in y (rows) and with all outlets on the 1st row.
        spacing = 300.0

        # top-border base level
        top_base_level = [
            NodeStatus.CORE,
            NodeStatus.CORE,
            NodeStatus.FIXED_VALUE,
            NodeStatus.CORE,
        ]

        grid = RasterGrid([2, 2], [spacing, spacing], top_base_level)

        flow_graph = FlowGraph(grid, [SingleFlowRouter()])

        area_exp = 0.5
        slope_exp = 1.0
        tolerance = 1e-3
        eroder = SPLEroder(flow_graph, k_coef, area_exp, slope_exp, tolerance)

        h = 1.0
        elevation = np.array([[0.0, 0.0], [h, h]], dtype="d")

        flow_graph.update_routes(elevation)
        drainage_area = flow_graph.accumulate(1.0)

        dt = 1  # use small time step (compare with explicit scheme)
        erosion = eroder.erode(elevation, drainage_area, dt)

        slope = h / spacing
        a = spacing**2
        k_coef = 1e-3
        err = dt * k_coef * a**area_exp * slope**slope_exp
        expected_erosion = np.array([[0.0, 0.0], [err, err]], dtype="d")

        np.testing.assert_allclose(erosion, expected_erosion, atol=1e-5)
        assert eroder.n_corr == 0


def _solve_diffusion_analytical(x, y, k_coef, t):
    fact = 4 * k_coef * t
    return np.exp(-(x * x + y * y) / fact) / (fact * np.pi)


def _compute_l2_norm(a1, a2):
    return 1.0 / a1.size * np.sum(a2**2 - a1**2)


class TestDiffusionADIEroder:
    def test_constructor_properties(self) -> None:
        grid = RasterGrid([2, 2], [1.0, 1.0], NodeStatus.FIXED_VALUE)

        eroder = DiffusionADIEroder(grid, 1e-3)

        np.testing.assert_equal(eroder.k_coef, np.full(grid.shape, 1e-3))

        k_coef_arr = np.full(grid.shape, 1e-2)
        eroder.k_coef = k_coef_arr
        np.testing.assert_equal(eroder.k_coef, k_coef_arr)

        k_coef = 1e-5
        eroder.k_coef = k_coef
        np.testing.assert_equal(eroder.k_coef, np.full(grid.shape, 1e-5))

        with pytest.raises(RuntimeError, match=".*shape mismatch"):
            eroder.k_coef = np.ones((4, 5))

    @pytest.mark.parametrize("k_coef_type", ["scalar", "array"])
    def test_erode_linear_diffusion(self, k_coef_type) -> None:
        x, y = np.meshgrid(np.linspace(-20, 20, 51), np.linspace(-20, 20, 101))
        dy = 0.4
        dx = 0.8

        grid = RasterGrid([101, 51], [dy, dx], NodeStatus.FIXED_VALUE)

        k_coef = 1e-3
        eroder = DiffusionADIEroder(grid, k_coef)

        dt = 1e3
        t0 = 2e3

        elevation_init = _solve_diffusion_analytical(x, y, k_coef, t0)
        erosion = np.empty_like(elevation_init)

        elevation_analytical = _solve_diffusion_analytical(x, y, k_coef, t0 + dt)

        if k_coef_type == "array":
            eroder.k_coef = np.full_like(x, k_coef)

        erosion = eroder.erode(elevation_init, dt)
        elevation_numerical = elevation_init - erosion

        l2_norm = _compute_l2_norm(elevation_analytical, elevation_numerical)

        assert l2_norm < 1e-9
