import numpy as np
import pytest

from fastscapelib.eroders import (
    erode_linear_diffusion_d,
    erode_linear_diffusion_var_d,
    erode_stream_power_d,
    erode_stream_power_var_d,
)
from fastscapelib.flow import (
    DummyFlowRouter,
    FlowGraph,
    MultipleFlowRouter,
    NoSinkResolver,
    SingleFlowRouter,
)
from fastscapelib.grid import (
    Node,
    NodeStatus,
    ProfileGrid,
    RasterBoundaryStatus,
    RasterGrid,
    RasterNode,
)


class TestErodeStreamPower:
    @pytest.mark.parametrize(
        "func,k",
        [(erode_stream_power_d, 1e-3), (erode_stream_power_var_d, np.full(4, 1e-3))],
    )
    def test_profile_grid(self, func, k):
        spacing = 300.0
        grid = ProfileGrid(4, 300, [NodeStatus.FIXED_VALUE_BOUNDARY] * 2, [])
        flow_graph = FlowGraph(grid, SingleFlowRouter(), NoSinkResolver())

        h = 1.0
        elevation = np.array([0.0, h, h, 0.0], dtype="d")
        graph_elevation = flow_graph.update_routes(elevation)

        drainage_area = flow_graph.accumulate(1.0)
        erosion = np.zeros_like(elevation)
        m_exp = 0.5
        n_exp = 1.0

        dt = 1.0  # use small time step (compare with explicit scheme)
        tolerance = 1e-3

        n_corr = func(
            erosion,
            elevation,
            drainage_area,
            flow_graph,
            k,
            m_exp,
            n_exp,
            dt,
            tolerance,
        )

        slope = h / spacing
        a = spacing
        k_coef = 1e-3
        err = dt * k_coef * a**m_exp * slope**n_exp
        expected_erosion = np.array([0.0, err, err, 0.0], dtype="d")

        np.testing.assert_allclose(erosion, expected_erosion, atol=1e-5)
        assert n_corr == 0

    @pytest.mark.parametrize(
        "func,k",
        [
            (erode_stream_power_d, 1e-3),
            (erode_stream_power_var_d, np.full((2, 2), 1e-3)),
        ],
    )
    def test_raster_grid(self, func, k):
        # Test on a tiny (2x2) 2-d square grid with a planar surface
        # tilted in y (rows) and with all outlets on the 1st row.
        spacing = 300.0
        grid = RasterGrid(
            [2, 2],
            [spacing, spacing],
            RasterBoundaryStatus(NodeStatus.FIXED_VALUE_BOUNDARY),
            [],
        )
        flow_graph = FlowGraph(grid, SingleFlowRouter(), NoSinkResolver())

        h = 1.0
        elevation = np.array([[0.0, 0.0], [h, h]], dtype="d")
        graph_elevation = flow_graph.update_routes(elevation)

        drainage_area = flow_graph.accumulate(1.0)
        erosion = np.zeros_like(elevation)
        m_exp = 0.5
        n_exp = 1.0

        dt = 1  # use small time step (compare with explicit scheme)
        tolerance = 1e-3

        n_corr = func(
            erosion,
            elevation,
            drainage_area,
            flow_graph,
            k,
            m_exp,
            n_exp,
            dt,
            tolerance,
        )

        slope = h / spacing
        a = spacing**2
        k_coef = 1e-3
        err = dt * k_coef * a**m_exp * slope**n_exp
        expected_erosion = np.array([[0.0, 0.0], [err, err]], dtype="d")

        np.testing.assert_allclose(erosion, expected_erosion, atol=1e-5)
        assert n_corr == 0


def _solve_diffusion_analytical(x, y, k_coef, t):
    fact = 4 * k_coef * t
    return np.exp(-(x * x + y * y) / fact) / (fact * np.pi)


def _compute_l2_norm(a1, a2):
    return 1.0 / a1.size * np.sum(a2**2 - a1**2)


@pytest.mark.parametrize("k_coef_type", ["constant", "variable"])
def test_erode_linear_diffusion(k_coef_type):
    x, y = np.meshgrid(np.linspace(-20, 20, 51), np.linspace(-20, 20, 101))
    dy = 0.4
    dx = 0.8

    k_coef = 1e-3
    dt = 1e3

    t0 = 2e3

    elevation_init = _solve_diffusion_analytical(x, y, k_coef, t0)
    erosion = np.empty_like(elevation_init)

    elevation_analytical = _solve_diffusion_analytical(x, y, k_coef, t0 + dt)

    if k_coef_type == "constant":
        func = erode_linear_diffusion_d
        k = k_coef
    elif k_coef_type == "variable":
        func = erode_linear_diffusion_var_d
        k = np.full_like(x, k_coef)

    func(erosion, elevation_init, k, dt, dx, dy)
    elevation_numerical = elevation_init - erosion

    l2_norm = _compute_l2_norm(elevation_analytical, elevation_numerical)

    assert l2_norm < 1e-9
