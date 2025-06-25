from typing import Any

import numpy as np
import pytest

from fastscapelib.eroders import (  # type: ignore[attr-defined]
    FlowKernelEroder,
    SPLEroder,
)
from fastscapelib.flow import FlowGraph, FlowGraphTraversalDir, SingleFlowRouter
from fastscapelib.grid import NodeStatus, RasterGrid

nb = pytest.importorskip("numba")


class BasicEroder(FlowKernelEroder):
    @staticmethod
    def param_spec():
        return dict()

    @staticmethod
    def input_spec():
        return dict(a=nb.float64, dt=nb.float64)

    @staticmethod
    def kernel_apply_dir():
        return FlowGraphTraversalDir.ANY

    @staticmethod
    def kernel_func(node):
        node.erosion = node.a * node.dt


def test_flow_kernel_eroder_basic() -> None:
    grid = RasterGrid([10, 10], [300.0, 300.0], NodeStatus.FIXED_VALUE)
    flow_graph = FlowGraph(grid, [SingleFlowRouter()])

    eroder = BasicEroder(flow_graph)

    actual = eroder.erode(a=3.0, dt=2.0)
    expected = np.ones(grid.shape) * 3 * 2
    np.testing.assert_array_equal(actual, expected)

    with pytest.raises(KeyError, match="inputs are missing"):
        eroder.erode(a=3.0)

    with pytest.raises(ValueError, match="invalid inputs"):
        eroder.erode(a=3.0, dt=2.0, not_an_input=0)


class InvalidEroder(FlowKernelEroder):
    @staticmethod
    def param_spec():
        return {}

    @staticmethod
    def input_spec():
        return {}

    @staticmethod
    def kernel_apply_dir():
        return FlowGraphTraversalDir.ANY

    @staticmethod
    def kernel_func(a, b):
        pass


def test_flow_kernel_eroder_invalid_kernel_func() -> None:
    grid = RasterGrid([10, 10], [300.0, 300.0], NodeStatus.FIXED_VALUE)
    flow_graph = FlowGraph(grid, [SingleFlowRouter()])

    with pytest.raises(TypeError, match="static method.*single argument"):
        InvalidEroder(flow_graph)


class SPLFlowKernelEroder(FlowKernelEroder):
    """Stream-Power Law implemented as a flow kernel eroder."""

    def __init__(
        self,
        flow_graph: FlowGraph,
        k_coef: float,
        area_exp: float,
        slope_exp: float,
        tolerance: float = 1e-3,
    ):
        super().__init__(flow_graph)

        self.kernel_data.bind(
            k_coef=k_coef, area_exp=area_exp, slope_exp=slope_exp, tolerance=tolerance
        )

    @staticmethod
    def param_spec():
        return {
            "k_coef": nb.float64,
            "area_exp": nb.float64,
            "slope_exp": nb.float64,
            "tolerance": nb.float64,
        }

    @staticmethod
    def input_spec():
        return {
            "elevation": nb.float64[::1],
            "drainage_area": nb.float64[::1],
            "dt": nb.float64,
        }

    @staticmethod
    def kernel_apply_dir() -> FlowGraphTraversalDir:
        return FlowGraphTraversalDir.BREADTH_UPSTREAM

    @staticmethod
    def kernel_func(node: Any):
        r_count = node.receivers.count
        if r_count == 1 and node.receivers.distance[0] == 0.0:
            return

        elevation_flooded = np.finfo(np.double).max

        for r in range(r_count):
            irec_elevation_next = (
                node.receivers.elevation[r] - node.receivers.erosion[r]
            )

            if irec_elevation_next < elevation_flooded:
                elevation_flooded = irec_elevation_next

        if node.elevation <= elevation_flooded:
            return

        eq_num = node.elevation
        eq_den = 1.0

        for r in range(r_count):
            irec_elevation = node.receivers.elevation[r]
            irec_elevation_next = irec_elevation - node.receivers.erosion[r]

            if irec_elevation > node.elevation:
                continue

            irec_weight = node.receivers.weight[r]
            irec_distance = node.receivers.distance[r]

            factor = (
                node.k_coef
                * node.dt
                * np.power(node.drainage_area * irec_weight, node.area_exp)
            )
            factor /= irec_distance
            eq_num += factor * irec_elevation_next
            eq_den += factor

        elevation_updated = eq_num / eq_den

        if elevation_updated < elevation_flooded:
            elevation_updated = elevation_flooded + np.finfo(np.double).tiny

        node.erosion = node.elevation - elevation_updated

    def erode(
        self, elevation: np.ndarray, drainage_area: np.ndarray, dt: float
    ) -> np.ndarray:
        return super().erode(elevation=elevation, drainage_area=drainage_area, dt=dt)


def test_spl_eroder_vs_kernel_eroder() -> None:
    grid = RasterGrid([10, 10], [300.0, 300.0], NodeStatus.FIXED_VALUE)
    flow_graph = FlowGraph(grid, [SingleFlowRouter()])

    eroder = SPLEroder(flow_graph, 1e-7, 0.4, 1, 1e-5)
    flow_kernel_eroder = SPLFlowKernelEroder(flow_graph, 1e-7, 0.4, 1, 1e-5)

    rng = np.random.Generator(np.random.PCG64(1234))
    init_elevation = rng.uniform(0, 5, size=grid.shape)
    elevation = init_elevation.copy()
    drainage_area = np.empty_like(init_elevation)
    uplift_rate = np.full_like(init_elevation, 1e-3)
    uplift_rate[[0, -1], :] = 0.0
    uplift_rate[:, [0, -1]] = 0.0

    # run a few time steps to test erosion array reset
    for dt in [1e4, 2e4]:
        uplift = dt * uplift_rate
        uplifted_elevation = elevation + uplift
        flow_graph.update_routes(uplifted_elevation)
        flow_graph.accumulate(drainage_area, 1.0)

        flow_kernel_spl_erosion = flow_kernel_eroder.erode(
            uplifted_elevation, drainage_area, dt
        )
        spl_erosion = eroder.erode(uplifted_elevation, drainage_area, dt)

        np.testing.assert_allclose(spl_erosion, flow_kernel_spl_erosion)

        elevation = uplifted_elevation - spl_erosion
