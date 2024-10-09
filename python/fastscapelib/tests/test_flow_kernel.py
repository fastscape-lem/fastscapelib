import numpy as np
import pytest

from fastscapelib.flow import (
    FlowGraph,
    FlowGraphTraversalDir,
    MultiFlowRouter,
    PFloodSinkResolver,
)
from fastscapelib.flow.numba import create_flow_kernel
from fastscapelib.grid import NodeStatus, RasterGrid

nb = pytest.importorskip("numba")


@pytest.fixture(scope="module")
def flow_graph():
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
    def kernel_func(_):
        pass

    yield kernel_func


@pytest.fixture(scope="module")
def compiled_kernel1(kernel_func1, flow_graph):
    kernel = create_flow_kernel(
        flow_graph,
        kernel_func1,
        spec=dict(
            a=nb.float64[::1],
        ),
        outputs=["a"],
        apply_dir=FlowGraphTraversalDir.ANY,
    )
    yield kernel


@pytest.fixture(scope="function")
def kernel1(compiled_kernel1):
    yield compiled_kernel1[0]


@pytest.fixture(scope="function")
def kernel1_data(compiled_kernel1):
    yield compiled_kernel1[1]
    # "un-bind" all kernel data
    compiled_kernel1[1]._bound_keys.clear()


@pytest.fixture(scope="module")
def compiled_kernel2(kernel_func1, flow_graph):
    kernel, data = create_flow_kernel(
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
        apply_dir=FlowGraphTraversalDir.ANY,
    )
    yield kernel, data


@pytest.fixture(scope="function")
def kernel2(compiled_kernel2):
    yield compiled_kernel2[0]


@pytest.fixture(scope="function")
def kernel2_data(compiled_kernel2):
    yield compiled_kernel2[1]
    # "un-bind" all kernel data
    compiled_kernel2[1]._bound_keys.clear()


@pytest.fixture(scope="module")
def compiled_kernel3(kernel_func2, flow_graph):
    kernel, data = create_flow_kernel(
        flow_graph,
        kernel_func2,
        spec=dict(
            a=nb.float64,
        ),
        apply_dir=FlowGraphTraversalDir.ANY,
        print_generated_code=True,
    )
    yield kernel, data


@pytest.fixture(scope="function")
def kernel3(compiled_kernel3):
    yield compiled_kernel3[0]


@pytest.fixture(scope="function")
def kernel3_data(compiled_kernel3):
    yield compiled_kernel3[1]
    # "un-bind" all kernel data
    compiled_kernel3[1]._bound_keys.clear()


class TestFlowKernelData:
    def test_bind_array(self, flow_graph, kernel1_data):
        # not bound yet
        assert kernel1_data.a is None

        # 1-d array case
        expected = np.zeros(flow_graph.size)
        kernel1_data.bind(a=expected)
        assert kernel1_data.a is expected

        # scalar (expand) case
        expected = np.ones(flow_graph.size)
        kernel1_data.bind(a=1.0)
        np.testing.assert_array_equal(kernel1_data.a, expected)

        # 2-d raster (flatten) case
        expected = 2 * np.ones(flow_graph.grid_shape)
        kernel1_data.bind(a=expected)
        np.testing.assert_array_equal(kernel1_data.a, expected.ravel())

        # invalid shape
        with pytest.raises(ValueError):
            kernel1_data.bind(a=np.zeros(flow_graph.size - 1))

    def test_bind_scalar(self, kernel2_data):
        # not bound yet
        assert kernel2_data.f64 is None
        kernel2_data.bind(f64=1.0)
        assert kernel2_data.f64 == 1.0

    def test_attr_like_access(self, kernel2_data):
        kernel2_data.bind(int64=10)
        assert kernel2_data.int64 == 10
        assert kernel2_data.jitclass_obj.int64 == 10

    def test_mapping_interface(self, kernel2_data):
        assert len(kernel2_data) == 11
        assert list(kernel2_data) == [
            "f64",
            "f32",
            "int32",
            "int64",
            "uint64",
            "f64_arr",
            "f32_arr",
            "int32_arr",
            "int64_arr",
            "uint64_arr",
            "a",
        ]

        # no data bound yet
        assert list(kernel2_data.values()) == [None] * 11

        kernel2_data.bind(f64=1.0)
        assert kernel2_data["f64"] == 1.0

    def test_multiple_bindings(self, flow_graph, kernel1_data):
        data = kernel1_data

        a1 = np.ones(flow_graph.size, dtype=np.float64)
        ref1_id = id(a1)
        data.bind(a=a1)
        assert ref1_id == id(data.a)
        np.testing.assert_almost_equal(
            data.a, np.ones(flow_graph.size, dtype=np.float64)
        )

        a2 = np.zeros(flow_graph.size, dtype=np.float64)
        ref2_id = id(a2)
        data.bind(a=a2)
        assert ref2_id == id(data.a)
        assert ref1_id != id(data.a)
        np.testing.assert_almost_equal(
            data.a, np.zeros(flow_graph.size, dtype=np.float64)
        )

    def test_check_data_bindings(self, flow_graph, kernel1, kernel1_data):
        with pytest.raises(ValueError):
            kernel1_data.check_bindings()

        with pytest.raises(ValueError):
            flow_graph.apply_kernel(kernel1, kernel1_data)

    def test_inline_bindings(self, flow_graph, kernel_func1):
        _, data = create_flow_kernel(
            flow_graph,
            kernel_func1,
            spec=dict(
                a=(nb.float32[::1], np.ones(flow_graph.size, dtype=np.float32) * 1.15),
            ),
            outputs=["a"],
            apply_dir=FlowGraphTraversalDir.ANY,
        )

        np.testing.assert_almost_equal(
            data.a, np.ones(flow_graph.size, dtype=np.float32) * 1.15
        )

    def test_ref_bindings(self, flow_graph, kernel_func1):
        a = np.ones(flow_graph.size, dtype=np.float32) * 1.15

        _, data = create_flow_kernel(
            flow_graph,
            kernel_func1,
            spec=dict(
                a=(nb.float32[::1], a),
            ),
            outputs=["a"],
            apply_dir=FlowGraphTraversalDir.ANY,
        )

        np.testing.assert_almost_equal(
            data.a, np.ones(flow_graph.size, dtype=np.float32) * 1.15
        )

        a *= 2.0
        np.testing.assert_almost_equal(
            data.a, np.ones(flow_graph.size, dtype=np.float32) * 2.3
        )


class TestFlowKernel:
    def test_input_assignment(self, flow_graph, kernel_func1):
        with pytest.raises(AttributeError):
            create_flow_kernel(
                flow_graph,
                kernel_func1,
                spec=dict(
                    a=nb.float64[::1],
                ),
                apply_dir=FlowGraphTraversalDir.ANY,
            )

    def test_output_assignment(self, flow_graph, kernel1, kernel1_data):
        kernel, data = kernel1, kernel1_data

        data.bind(a=np.zeros(flow_graph.size))
        flow_graph.apply_kernel(kernel, data)

        np.testing.assert_almost_equal(
            data.a, np.ones(flow_graph.size, dtype=np.float64) * 10.0
        )

    def test_scalar_output(self, flow_graph, kernel_func1):
        with pytest.raises(TypeError):
            create_flow_kernel(
                flow_graph,
                kernel_func1,
                spec=dict(
                    a=(nb.float64, 10.0),
                ),
                outputs=["a"],
                apply_dir=FlowGraphTraversalDir.ANY,
            )

    def test_invalid_output(self, flow_graph, kernel_func1):
        with pytest.raises(KeyError):
            create_flow_kernel(
                flow_graph,
                kernel_func1,
                spec=dict(
                    a=(nb.float64, 10.0),
                ),
                outputs=["b"],
                apply_dir=FlowGraphTraversalDir.ANY,
            )

    def test_multiple_types(self, kernel2, kernel2_data):
        kernel, data = kernel2, kernel2_data

        assert data.jitclass_obj._numba_type_.struct["f64"] == nb.float64
        assert data.jitclass_obj._numba_type_.struct["f32"] == nb.float32
        assert data.jitclass_obj._numba_type_.struct["int32"] == nb.int32
        assert data.jitclass_obj._numba_type_.struct["int64"] == nb.int64
        assert data.jitclass_obj._numba_type_.struct["uint64"] == nb.uint64

        assert data.jitclass_obj._numba_type_.struct["f64_arr"] == nb.float64[::1]
        assert data.jitclass_obj._numba_type_.struct["f32_arr"] == nb.float32[::1]
        assert data.jitclass_obj._numba_type_.struct["int32_arr"] == nb.int32[::1]
        assert data.jitclass_obj._numba_type_.struct["int64_arr"] == nb.int64[::1]
        assert data.jitclass_obj._numba_type_.struct["uint64_arr"] == nb.uint64[::1]

        assert type(data.jitclass_obj.f64_arr) == np.ndarray
        assert data.jitclass_obj.f64_arr.dtype == np.float64
        assert type(data.jitclass_obj.f32_arr) == np.ndarray
        assert data.jitclass_obj.f32_arr.dtype == np.float32
        assert type(data.jitclass_obj.int32_arr) == np.ndarray
        assert data.jitclass_obj.int32_arr.dtype == np.int32
        assert type(data.jitclass_obj.int64_arr) == np.ndarray
        assert data.jitclass_obj.int64_arr.dtype == np.int64
        assert type(data.jitclass_obj.uint64_arr) == np.ndarray
        assert data.jitclass_obj.uint64_arr.dtype == np.uint64

        node_data = kernel.node_data_create()
        node_data_struct = node_data.__class__._numba_type_.class_type.struct
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
        node_data = kernel1.node_data_create()
        node_data_struct = node_data.__class__._numba_type_.class_type.struct
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
        assert node_data.__class__.__name__ == "FlowKernelNodeData"
        assert hasattr(node_data, "_numba_type_")
        assert hasattr(node_data, "a")

    def test_node_data_init(self, kernel1, kernel3, kernel3_data):
        assert kernel1.node_data_init is None
        assert kernel3.node_data_init is not None

        node_data = kernel3.node_data_create()
        kernel3_data.bind(a=12.0)
        kernel3.node_data_init(node_data, kernel3_data.jitclass_obj)
        assert node_data.a == 12.0

        kernel3_data.bind(a=1.0)
        kernel3.node_data_init(node_data, kernel3_data.jitclass_obj)
        assert node_data.a == 1.0

    def test_max_receivers(self, flow_graph, kernel1, kernel1_data, kernel_func1):
        rng = np.random.Generator(np.random.PCG64(1234))
        init_elevation = rng.uniform(0, 5, size=flow_graph.grid_shape)

        flow_graph.update_routes(init_elevation)
        assert np.any(flow_graph.impl().receivers_count > 1)

        kernel1_data.bind(a=np.ones(flow_graph.size, dtype=np.float64))
        flow_graph.apply_kernel(kernel1, kernel1_data)

        node_data = kernel1.node_data_create()
        for i in range(flow_graph.size):
            kernel1.node_data_getter(i, kernel1_data.jitclass_obj, node_data)
            assert node_data.receivers.count == flow_graph.impl().receivers_count[i]

        kernel, data = create_flow_kernel(
            flow_graph,
            kernel_func1,
            spec=dict(
                a=nb.float64[::1],
            ),
            outputs=["a"],
            apply_dir=FlowGraphTraversalDir.ANY,
            max_receivers=1,
        )
        data.bind(a=np.ones(flow_graph.size, dtype=np.float64))
        with pytest.raises(ValueError):
            flow_graph.apply_kernel(kernel, data)
