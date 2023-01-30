import numpy as np
import numpy.testing as npt
import pytest
from fastscapelib.flow import (
    FlowDirection,
    FlowGraph,
    FlowSnapshot,
    MSTSinkResolver,
    PFloodSinkResolver,
    SingleFlowRouter,
)
from fastscapelib.grid import NodeStatus, RasterBoundaryStatus, RasterGrid


@pytest.fixture
def grid():
    return RasterGrid(
        [10, 20],
        [1.0, 2.0],
        RasterBoundaryStatus(NodeStatus.FIXED_VALUE_BOUNDARY),
        [],
    )


def test_flow_snapshot_class_attr():
    assert FlowSnapshot.graph_updated is False
    assert FlowSnapshot.elevation_updated is False
    assert FlowSnapshot.in_flowdir == FlowDirection.UNDEFINED
    assert FlowSnapshot.out_flowdir == FlowDirection.UNDEFINED


def test_flow_snapshot_error(grid):
    with pytest.raises(ValueError, match="no flow routing operator defined before"):
        FlowGraph(grid, [FlowSnapshot("test"), SingleFlowRouter()])


def test_flow_snapshot_items(grid):
    graph = FlowGraph(
        grid,
        [
            PFloodSinkResolver(),
            FlowSnapshot("a", save_graph=False, save_elevation=True),
            SingleFlowRouter(),
            FlowSnapshot("b", save_graph=True, save_elevation=False),
        ],
    )

    assert graph.graph_snapshot_keys == ["b"]
    assert graph.elevation_snapshot_keys == ["a"]

    assert isinstance(graph.graph_snapshot("b"), FlowGraph)
    assert isinstance(graph.elevation_snapshot("a"), np.ndarray)


def test_snapshot_graph(grid):
    graph = FlowGraph(
        grid,
        [SingleFlowRouter(), FlowSnapshot("a"), MSTSinkResolver(), FlowSnapshot("b")],
    )

    elevation = np.random.uniform(size=grid.shape)
    graph.update_routes(elevation)
    snapshot_a = graph.graph_snapshot("a")
    snapshot_b = graph.graph_snapshot("b")

    attrs = [
        "receivers_count",
        "receivers",
        "receivers_distance",
        "receivers_weight",
        "donors_count",
        "dfs_indices",
    ]

    for name in attrs:
        arr = getattr(graph.impl(), name)
        snapshot_a_arr = getattr(snapshot_a.impl(), name)
        snapshot_b_arr = getattr(snapshot_b.impl(), name)

        assert arr is not snapshot_a_arr

        npt.assert_equal(arr, snapshot_b_arr)
        assert arr is not snapshot_b_arr

    npt.assert_equal(snapshot_b.impl().donors[:, 0], graph.impl().donors[:, 0])
    assert not np.all(snapshot_a.impl().receivers == snapshot_b.impl().receivers)

    # high-level interface test
    actual = graph.accumulate(1.0)
    expected = snapshot_b.accumulate(1.0)
    npt.assert_allclose(actual, expected)

    # test snapshot graphs are read-only
    with pytest.raises(RuntimeError, match=".*(read-only)"):
        snapshot_a.update_routes(elevation)


def test_snapshot_elevation(grid):
    graph = FlowGraph(
        grid,
        [
            FlowSnapshot("a", save_graph=False, save_elevation=True),
            PFloodSinkResolver(),
            FlowSnapshot("b", save_graph=False, save_elevation=True),
            SingleFlowRouter(),
        ],
    )

    elevation = np.random.uniform(size=grid.shape)
    filled_elevation = graph.update_routes(elevation)
    snapshot_a = graph.elevation_snapshot("a")
    snapshot_b = graph.elevation_snapshot("b")

    npt.assert_equal(elevation, snapshot_a)
    assert elevation is not snapshot_a
    npt.assert_equal(filled_elevation, snapshot_b)
    assert filled_elevation is not snapshot_b
