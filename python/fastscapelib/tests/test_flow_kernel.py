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
def grid():
    raster_grid = RasterGrid([5, 5], [2.2, 2.4], NodeStatus.FIXED_VALUE)
    yield raster_grid


@pytest.fixture(scope="module")
def flow_graph(grid):
    flow_graph = FlowGraph(grid, [PFloodSinkResolver(), MultiFlowRouter()])
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
    kernel, data = create_flow_kernel(
        flow_graph,
        kernel_func1,
        spec=dict(
            a=nb.float64[::1],
        ),
        outputs=["a"],
        apply_dir=FlowGraphTraversalDir.ANY,
    )
    yield kernel, data


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
            b=nb.float64[::1],
        ),
        outputs=["b"],
        apply_dir=FlowGraphTraversalDir.ANY,
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

    def test_var_names(self, kernel2_data):
        assert kernel2_data.var_names == tuple(kernel2_data)
        assert kernel2_data.grid_var_names == (
            "f64_arr",
            "f32_arr",
            "int32_arr",
            "int64_arr",
            "uint64_arr",
            "a",
        )

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
    def test_no_output_error(self, flow_graph, kernel_func1):
        with pytest.raises(ValueError, match="no output variable"):
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

    def test_scalar_output_error(self, flow_graph, kernel_func1):
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

    def test_output_not_in_spec_error(self, flow_graph, kernel_func1):
        with pytest.raises(ValueError, match=".*output variables.*not defined in spec"):
            create_flow_kernel(
                flow_graph,
                kernel_func1,
                spec=dict(
                    a=(nb.float64, 10.0),
                ),
                outputs=["b"],
                apply_dir=FlowGraphTraversalDir.ANY,
            )

    def test_kernel_properties(self, kernel1):
        assert kernel1.n_threads == 1
        assert kernel1.apply_dir == FlowGraphTraversalDir.ANY
        assert kernel1.outputs == ("a",)

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

        node_data_donors_struct = node_data_struct["donors"].struct
        assert len(node_data_donors_struct) == 13
        assert node_data_donors_struct["a"] == nb.float64[::1]
        assert node_data_donors_struct["f64_arr"] == nb.float64[::1]
        assert node_data_donors_struct["f32_arr"] == nb.float32[::1]
        assert node_data_donors_struct["int32_arr"] == nb.int32[::1]
        assert node_data_donors_struct["int64_arr"] == nb.int64[::1]
        assert node_data_donors_struct["uint64_arr"] == nb.uint64[::1]
        assert node_data_donors_struct["count"] == nb.uint64

    def test_node_data(self, kernel1):
        node_data = kernel1.node_data_create()
        node_data_struct = node_data.__class__._numba_type_.class_type.struct
        assert len(node_data_struct) == 3
        assert "receivers" in node_data_struct
        assert "donors" in node_data_struct
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

        node_data_donors_struct = node_data_struct["donors"].struct
        assert len(node_data_donors_struct) == 3
        assert node_data_donors_struct["a"] == nb.float64[::1]
        assert node_data_donors_struct["_a"] == nb.float64[::1]
        assert node_data_donors_struct["count"] == nb.uint64

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


def test_reserved_spec(flow_graph):
    def kernel_func(_): ...

    with pytest.raises(ValueError, match="reserved variable names defined in spec"):
        create_flow_kernel(
            flow_graph,
            kernel_func,
            spec=dict(receivers_count=nb.int64[::1]),
            outputs=["receivers_count"],
        )


def test_simple_recounter_flow_kernel(grid, flow_graph):
    # recount flow receivers and donors at each node in a
    # flow kernel function by reading values at the receivers and donors
    # and updating output values at the current node.
    def recounter_kernel_func(node):
        nb_receivers = 0
        nb_donors = 0

        for i in range(node.receivers.count):
            nb_receivers += node.receivers.one[i]
        for i in range(node.donors.count):
            nb_donors += node.donors.one[i]

        node.rec_count = nb_receivers
        node.don_count = nb_donors

    kernel, data = create_flow_kernel(
        flow_graph,
        recounter_kernel_func,
        spec=dict(
            one=nb.int64[::1],
            rec_count=nb.int64[::1],
            don_count=nb.int64[::1],
        ),
        outputs=["rec_count", "don_count"],
        # FIXME: output values not updated with auto_resize=True (dynamic resizing)
        auto_resize=False,
        n_threads=1,
        apply_dir=FlowGraphTraversalDir.ANY,
    )

    elevation = np.random.uniform(size=grid.shape)
    flow_graph.update_routes(elevation)

    one = np.ones(grid.size, dtype=np.int64)
    rec_count = np.zeros(grid.size, dtype=np.int64)
    don_count = np.zeros(grid.size, dtype=np.int64)
    data.bind(one=one, rec_count=rec_count, don_count=don_count)
    flow_graph.apply_kernel(kernel, data)

    np.testing.assert_array_equal(rec_count, flow_graph.impl().receivers_count)
    np.testing.assert_array_equal(don_count, flow_graph.impl().donors_count)


def test_simple_recounter_flow_kernel2(grid, flow_graph):
    # recount flow receivers and donors at each node in a
    # flow kernel function by directly updating output values at
    # the receivers and donors of the current node.
    def recounter_kernel_func2(node):
        for i in range(node.donors.count):
            node.donors.rec_count[i] += 1

        if node.receivers.count == 1 and node.receivers.distance[0] == 0.0:
            # base level node
            node.rec_count = 1
        else:
            for i in range(node.receivers.count):
                node.receivers.don_count[i] += 1

    kernel, data = create_flow_kernel(
        flow_graph,
        recounter_kernel_func2,
        spec=dict(
            rec_count=nb.int64[::1],
            don_count=nb.int64[::1],
        ),
        outputs=["rec_count", "don_count"],
        # FIXME: output values not updated with auto_resize=True (dynamic resizing)
        auto_resize=False,
        n_threads=1,
        apply_dir=FlowGraphTraversalDir.ANY,
    )

    elevation = np.random.uniform(size=grid.shape)
    flow_graph.update_routes(elevation)

    rec_count = np.zeros(grid.size, dtype=np.int64)
    don_count = np.zeros(grid.size, dtype=np.int64)
    data.bind(rec_count=rec_count, don_count=don_count)
    flow_graph.apply_kernel(kernel, data)

    np.testing.assert_array_equal(rec_count, flow_graph.impl().receivers_count)
    np.testing.assert_array_equal(don_count, flow_graph.impl().donors_count)


@pytest.mark.parametrize(
    "apply_dir",
    [FlowGraphTraversalDir.DEPTH_DOWNSTREAM, FlowGraphTraversalDir.BREADTH_DOWNSTREAM],
)
def test_simple_accumulate_flow_kernel(grid, flow_graph, apply_dir):
    def accumulate_kernel_func(node):
        r_count = node.receivers.count
        if r_count == 1 and node.receivers.distance[0] == 0.0:
            # base level node
            return
        for r in range(r_count):
            irec_weight = node.receivers.weight[r]
            node.receivers.drainage_area[r] += node.drainage_area * irec_weight

    kernel, data = create_flow_kernel(
        flow_graph,
        accumulate_kernel_func,
        spec=dict(
            drainage_area=nb.float64[::1],
        ),
        outputs=["drainage_area"],
        # FIXME: output values not updated with auto_resize=True (dynamic resizing)
        auto_resize=False,
        n_threads=1,
        apply_dir=apply_dir,
    )

    elevation = np.random.uniform(size=grid.shape)
    flow_graph.update_routes(elevation)

    drainage_area_actual = grid.nodes_areas()
    data.bind(drainage_area=drainage_area_actual.ravel())
    flow_graph.apply_kernel(kernel, data)

    drainage_area_expected = flow_graph.accumulate(1.0)

    np.testing.assert_allclose(drainage_area_actual, drainage_area_expected)


@pytest.mark.parametrize(
    "receiver_or_donor,data_access",
    [
        ("receiver", True),
        ("receiver", False),
        ("donor", True),
        ("donor", False),
    ],
)
def test_data_access_at_receivers_and_donors(
    flow_graph, receiver_or_donor, data_access
):
    def kernel_func(_): ...

    kernel, _ = create_flow_kernel(
        flow_graph,
        kernel_func,
        spec=dict(a=nb.float64[::1]),
        outputs=["a"],
        get_data_at_receivers=data_access,
        set_data_at_receivers=data_access,
        get_data_at_donors=data_access,
        set_data_at_donors=data_access,
    )

    # receivers and donors content node data always generated to prevent
    # segmentation faults
    line_content = f"self.{receiver_or_donor}s.a = self.{receiver_or_donor}s._a[:]"
    assert line_content in kernel.generated_code["node_data_jitclass_init"]

    line_view = f"({receiver_or_donor}s.a, {receiver_or_donor}s._a)"
    line_content = f"{receiver_or_donor}s._a[i] = data.a[{receiver_or_donor}_idx]"
    if data_access:
        assert line_view in kernel.generated_code["node_data_getter"]
        assert line_content in kernel.generated_code["node_data_getter"]
    else:
        assert line_view not in kernel.generated_code["node_data_getter"]
        assert line_content not in kernel.generated_code["node_data_getter"]

    line_content = f"data.a[{receiver_or_donor}_idx] = {receiver_or_donor}s.a[i]"
    if data_access:
        assert line_content in kernel.generated_code["node_data_setter"]
    else:
        assert line_content not in kernel.generated_code["node_data_setter"]
