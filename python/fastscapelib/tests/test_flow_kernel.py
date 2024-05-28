from textwrap import dedent

import numba as nb
import numpy as np
import numpy.testing as npt
import pytest

from fastscapelib.flow import (
    FlowGraph,
    FlowSnapshot,
    KernelApplicationOrder,
    MultiFlowRouter,
    NumbaFlowKernel,
    PFloodSinkResolver,
    SingleFlowRouter,
)
from fastscapelib.grid import NodeStatus, ProfileGrid, RasterGrid


@pytest.fixture(scope="module")
def flow_graph() -> None:
    raster_grid = RasterGrid([5, 5], [2.2, 2.4], NodeStatus.FIXED_VALUE)
    flow_graph = FlowGraph(raster_grid, [PFloodSinkResolver(), MultiFlowRouter()])
    yield flow_graph


@pytest.fixture(scope="module")
def kernel_func1():
    def kernel_func(node):
        node.a = 10.0

    yield kernel_func


@pytest.fixture(scope="module")
def kernel_func2():
    def kernel_func(node):
        pass

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


@pytest.fixture(scope="module")
def compiled_kernel2(kernel_func1, flow_graph):
    kernel = NumbaFlowKernel(
        flow_graph,
        kernel_func1,
        spec=dict(
            f64=nb.float64,
            f32=nb.float32,
            int32=nb.int32,
            int64=nb.int64,
            uint64=nb.uint64,
            f64_arr=nb.float64[::1],
            f32_arr=nb.float32[::1],
            int32_arr=nb.int32[::1],
            int64_arr=nb.int64[::1],
            uint64_arr=nb.uint64[::1],
            a=nb.float64[::1],
        ),
        outputs=["a"],
        application_order=KernelApplicationOrder.ANY,
    )
    yield kernel


@pytest.fixture(scope="function")
def kernel2(compiled_kernel2):
    yield compiled_kernel2

    compiled_kernel2._bound_data.clear()


@pytest.fixture(scope="module")
def compiled_kernel3(kernel_func2, flow_graph):
    kernel = NumbaFlowKernel(
        flow_graph,
        kernel_func2,
        spec=dict(
            a=nb.float64,
        ),
        application_order=KernelApplicationOrder.ANY,
        print_generated_code=True,
    )
    yield kernel


@pytest.fixture(scope="function")
def kernel3(compiled_kernel3):
    yield compiled_kernel3

    compiled_kernel3._bound_data.clear()


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

        np.testing.assert_almost_equal(
            kernel1._data.a, np.ones(flow_graph.size, dtype=np.float64) * 10.0
        )

    def test_multiple_bindings(self, flow_graph, kernel1):

        a1 = np.ones(flow_graph.size, dtype=np.float64)
        ref1_id = id(a1)
        kernel1.bind_data(a=a1)
        assert ref1_id == id(kernel1._data.a)
        np.testing.assert_almost_equal(
            kernel1._data.a, np.ones(flow_graph.size, dtype=np.float64)
        )

        a2 = np.zeros(flow_graph.size, dtype=np.float64)
        ref2_id = id(a2)
        kernel1.bind_data(a=a2)
        assert ref2_id == id(kernel1._data.a)
        assert ref1_id != id(kernel1._data.a)
        np.testing.assert_almost_equal(
            kernel1._data.a, np.zeros(flow_graph.size, dtype=np.float64)
        )

    def test_check_data_bindings(self, flow_graph, kernel1):
        with pytest.raises(ValueError):
            kernel1.check_data_bindings()

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

    def test_scalar_output(self, flow_graph, kernel_func1):
        with pytest.raises(TypeError):
            NumbaFlowKernel(
                flow_graph,
                kernel_func1,
                spec=dict(
                    a=10.0,
                ),
                outputs=["a"],
                application_order=KernelApplicationOrder.ANY,
            )

    def test_invalid_output(self, flow_graph, kernel_func1):
        with pytest.raises(KeyError):
            NumbaFlowKernel(
                flow_graph,
                kernel_func1,
                spec=dict(
                    a=10.0,
                ),
                outputs=["b"],
                application_order=KernelApplicationOrder.ANY,
            )

    def test_invalid_grid_data_shape(self, flow_graph, kernel_func1):
        with pytest.raises(AttributeError):
            kernel1.bind_data(a=np.zeros(flow_graph.size - 1))

    def test_multiple_types(self, kernel2):

        assert kernel2._data._numba_type_.struct["f64"] == nb.float64
        assert kernel2._data._numba_type_.struct["f32"] == nb.float32
        assert kernel2._data._numba_type_.struct["int32"] == nb.int32
        assert kernel2._data._numba_type_.struct["int64"] == nb.int64
        assert kernel2._data._numba_type_.struct["uint64"] == nb.uint64

        assert kernel2._data._numba_type_.struct["f64_arr"] == nb.float64[::1]
        assert kernel2._data._numba_type_.struct["f32_arr"] == nb.float32[::1]
        assert kernel2._data._numba_type_.struct["int32_arr"] == nb.int32[::1]
        assert kernel2._data._numba_type_.struct["int64_arr"] == nb.int64[::1]
        assert kernel2._data._numba_type_.struct["uint64_arr"] == nb.uint64[::1]

        assert type(kernel2._data.f64_arr) == np.ndarray
        assert kernel2._data.f64_arr.dtype == np.float64
        assert type(kernel2._data.f32_arr) == np.ndarray
        assert kernel2._data.f32_arr.dtype == np.float32
        assert type(kernel2._data.int32_arr) == np.ndarray
        assert kernel2._data.int32_arr.dtype == np.int32
        assert type(kernel2._data.int64_arr) == np.ndarray
        assert kernel2._data.int64_arr.dtype == np.int64
        assert type(kernel2._data.uint64_arr) == np.ndarray
        assert kernel2._data.uint64_arr.dtype == np.uint64

        node_data_struct = kernel2._node_data_jitclass.class_type.struct
        assert node_data_struct["f64_arr"] == nb.float64
        assert node_data_struct["f32_arr"] == nb.float32
        assert node_data_struct["int32_arr"] == nb.int32
        assert node_data_struct["int64_arr"] == nb.int64
        assert node_data_struct["uint64_arr"] == nb.uint64

        node_data_receivers_struct = node_data_struct["receivers"].struct
        assert len(node_data_receivers_struct) == 17
        assert node_data_receivers_struct["a"] == nb.float64[::1]
        assert node_data_receivers_struct["f64_arr"] == nb.float64[::1]
        assert node_data_receivers_struct["f32_arr"] == nb.float32[::1]
        assert node_data_receivers_struct["int32_arr"] == nb.int32[::1]
        assert node_data_receivers_struct["int64_arr"] == nb.int64[::1]
        assert node_data_receivers_struct["uint64_arr"] == nb.uint64[::1]
        assert node_data_receivers_struct["distance"] == nb.float64[::1]
        assert node_data_receivers_struct["weight"] == nb.float64[::1]
        assert node_data_receivers_struct["count"] == nb.uint64

    def test_node_data(self, kernel1):

        node_data_struct = kernel1._node_data_jitclass.class_type.struct
        assert len(node_data_struct) == 2
        assert "receivers" in node_data_struct
        assert node_data_struct["a"] == nb.float64

        node_data_receivers_struct = node_data_struct["receivers"].struct
        assert len(node_data_receivers_struct) == 7
        assert node_data_receivers_struct["a"] == nb.float64[::1]
        assert node_data_receivers_struct["_a"] == nb.float64[::1]
        assert node_data_receivers_struct["distance"] == nb.float64[::1]
        assert node_data_receivers_struct["_distance"] == nb.float64[::1]
        assert node_data_receivers_struct["weight"] == nb.float64[::1]
        assert node_data_receivers_struct["_weight"] == nb.float64[::1]
        assert node_data_receivers_struct["count"] == nb.uint64

    def test_node_data_create(self, kernel1):

        node_data = kernel1.node_data_create()
        assert (
            node_data._numba_type_
            is kernel1._node_data_jitclass.class_type.instance_type
        )

    def test_node_data_init(self, kernel1, kernel3):

        assert kernel1.node_data_init is None
        assert kernel3.node_data_init is not None

        node_data = kernel3.node_data_create()
        kernel3._data.a = 12.0
        kernel3.node_data_init(node_data, kernel3._data)
        assert node_data.a == 12.0

        kernel3._data.a = 1.0
        kernel3.node_data_init(node_data, kernel3._data)
        assert node_data.a == 1.0

    def test_max_receivers(self, flow_graph, kernel1, kernel_func1):
        rng = np.random.Generator(np.random.PCG64(1234))
        init_elevation = rng.uniform(0, 5, size=flow_graph.grid_shape)

        flow_graph.update_routes(init_elevation)
        assert np.any(flow_graph.impl().receivers_count > 1)

        kernel1.bind_data(a=np.ones(flow_graph.size, dtype=np.float64))
        flow_graph.apply_kernel(kernel1, kernel1.data)

        node_data = kernel1.node_data_create()
        for i in range(flow_graph.size):
            kernel1.node_data_getter(i, kernel1._data, node_data)
            assert node_data.receivers.count == flow_graph.impl().receivers_count[i]

        kernel = NumbaFlowKernel(
            flow_graph,
            kernel_func1,
            spec=dict(
                a=nb.float64[::1],
            ),
            outputs=["a"],
            application_order=KernelApplicationOrder.ANY,
            max_receivers=1,
        )
        kernel.bind_data(a=np.ones(flow_graph.size, dtype=np.float64))
        with pytest.raises(RuntimeError):
            flow_graph.apply_kernel(kernel, kernel.data)
