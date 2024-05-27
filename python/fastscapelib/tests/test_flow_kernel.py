from textwrap import dedent

import numpy as np
import numpy.testing as npt
import pytest

from fastscapelib.flow import (
    FlowGraph,
    FlowSnapshot,
    MultiFlowRouter,
    PFloodSinkResolver,
    SingleFlowRouter,
    NumbaFlowKernel,
    KernelApplicationOrder,
)
from fastscapelib.grid import NodeStatus, ProfileGrid, RasterGrid
import numba as nb


@pytest.fixture(scope="module")
def flow_graph() -> None:
    raster_grid = RasterGrid([5, 10], [2.2, 2.4], NodeStatus.FIXED_VALUE)
    flow_graph = FlowGraph(raster_grid, [PFloodSinkResolver(), SingleFlowRouter()])
    yield flow_graph


@pytest.fixture(scope="module")
def kernel_func1():
    def kernel_func(node):
        node.a = 10.0

    yield kernel_func


@pytest.fixture(scope="module")
def compiled_kernel1(kernel_func1, flow_graph):
    kernel = NumbaFlowKernel(
        flow_graph,
        kernel_func1,
        spec=dict(
            a=nb.float64[::1],
        ),
        outputs=["a"],
        application_order=KernelApplicationOrder.ANY,
    )
    yield kernel


@pytest.fixture(scope="function")
def kernel1(compiled_kernel1):
    yield compiled_kernel1

    compiled_kernel1._bound_data.clear()


class TestFlowKernel:

    def test_input_assignment(self, flow_graph, kernel_func1):
        with pytest.raises(AttributeError):
            NumbaFlowKernel(
                flow_graph,
                kernel_func1,
                spec=dict(
                    a=nb.float64[::1],
                ),
                application_order=KernelApplicationOrder.ANY,
            )

    def test_output_assignment(self, flow_graph, kernel1):
        assert "a" in kernel1._grid_data_ty
        assert len(kernel1._scalar_data_ty) == 0

        kernel1.bind_data(a=np.zeros(flow_graph.size))
        flow_graph.apply_kernel(kernel1, kernel1.data)

    def test_bindings_check(self, flow_graph, kernel1):
        assert "a" in kernel1._grid_data_ty
        assert len(kernel1._scalar_data_ty) == 0

        with pytest.raises(ValueError):
            flow_graph.apply_kernel(kernel1, kernel1.data)

    def test_inline_bindings(self, flow_graph, kernel_func1):

        kernel = NumbaFlowKernel(
            flow_graph,
            kernel_func1,
            spec=dict(
                a=np.ones(flow_graph.size, dtype=np.float32) * 1.15,
            ),
            outputs=["a"],
            application_order=KernelApplicationOrder.ANY,
        )

        np.testing.assert_almost_equal(
            kernel._data.a, np.ones(flow_graph.size, dtype=np.float32) * 1.15
        )

    def test_ref_bindings(self, flow_graph, kernel_func1):
        a = np.ones(flow_graph.size, dtype=np.float32) * 1.15

        kernel = NumbaFlowKernel(
            flow_graph,
            kernel_func1,
            spec=dict(
                a=a,
            ),
            outputs=["a"],
            application_order=KernelApplicationOrder.ANY,
        )

        np.testing.assert_almost_equal(
            kernel._data.a, np.ones(flow_graph.size, dtype=np.float32) * 1.15
        )

        a *= 2.0
        np.testing.assert_almost_equal(
            kernel._data.a, np.ones(flow_graph.size, dtype=np.float32) * 2.3
        )
