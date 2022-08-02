import numpy as np
import pytest

from fastscapelib.flow import (
    FlowGraph,
    MSTMethod,
    MSTRouteMethod,
    MSTSinkResolver,
    NoSinkResolver,
    PFloodSinkResolver,
    SingleFlowRouter,
)
from fastscapelib.grid import NodeStatus, RasterBoundaryStatus, RasterGrid


@pytest.fixture
def grid():
    bs = RasterBoundaryStatus(
        [
            NodeStatus.CORE,
            NodeStatus.CORE,
            NodeStatus.CORE,
            NodeStatus.FIXED_VALUE_BOUNDARY,
        ]
    )
    grid = RasterGrid([5, 5], [1, 1], bs, [])

    yield grid


@pytest.fixture
def elevation():
    # planar surface tilted along the y-axis
    # + one closed depression with pit at row/col index (1, 3)
    #   and pass between (2, 2) and (3, 1)
    elevation = np.array(
        [
            [1.00, 1.00, 1.00, 1.00, 1.00],
            [0.70, 0.70, 0.70, 0.10, 0.70],
            [0.50, 0.50, 0.15, 0.50, 0.50],
            [0.20, 0.19, 0.20, 0.20, 0.20],
            [0.00, 0.00, 0.00, 0.00, 0.00],
        ]
    )

    yield elevation


def test_no_sink_resolver(grid, elevation):
    graph = FlowGraph(grid, SingleFlowRouter(), NoSinkResolver())

    new_elevation = graph.update_routes(elevation)

    # unchanged elevation
    np.testing.assert_array_equal(new_elevation, elevation)

    # flow trapped in pit
    pit_idx_flat = 8
    assert graph.impl().receivers[pit_idx_flat, 0] == pit_idx_flat


def test_pflood_sink_resolver(grid, elevation):
    graph = FlowGraph(grid, SingleFlowRouter(), PFloodSinkResolver())

    new_elevation = graph.update_routes(elevation)

    # filled elevation (tiny slope)
    mask_filled = np.zeros_like(elevation, dtype=bool)
    mask_filled[[1, 2], [3, 2]] = True
    np.testing.assert_array_equal(new_elevation[~mask_filled], elevation[~mask_filled])
    assert new_elevation[2, 2] > new_elevation[3, 1]
    assert new_elevation[1, 3] > new_elevation[2, 2]

    # resolved flow
    assert graph.impl().receivers[8, 0] == 12
    assert graph.impl().receivers_distance[8, 0] == np.sqrt(2)
    assert graph.impl().receivers[12, 0] == 16
    assert graph.impl().receivers_distance[12, 0] == np.sqrt(2)


class TestMSTSinkResolver:
    def test_constructor(self):
        resolver = MSTSinkResolver()
        assert resolver.basin_method == MSTMethod.KRUSKAL
        assert resolver.route_method == MSTRouteMethod.CARVE

        # read-write attributes
        resolver.basin_method = MSTMethod.BORUVKA
        assert resolver.basin_method == MSTMethod.BORUVKA

        resolver.route_method = MSTRouteMethod.BASIC
        assert resolver.route_method == MSTRouteMethod.BASIC

        # constructor with arguments
        resolver2 = MSTSinkResolver(MSTMethod.BORUVKA, MSTRouteMethod.BASIC)
        assert resolver2.basin_method == MSTMethod.BORUVKA
        assert resolver2.route_method == MSTRouteMethod.BASIC

    @pytest.mark.parametrize("mst_method", [MSTMethod.KRUSKAL, MSTMethod.BORUVKA])
    def test_resolve_basic(self, grid, elevation, mst_method):

        resolver = MSTSinkResolver(mst_method, MSTRouteMethod.BASIC)
        graph = FlowGraph(grid, SingleFlowRouter(), resolver)

        new_elevation = graph.update_routes(elevation)

        # filled elevation (tiny slope)
        mask_filled = np.zeros_like(elevation, dtype=bool)
        mask_filled[[1, 2], [3, 2]] = True
        np.testing.assert_array_equal(
            new_elevation[~mask_filled], elevation[~mask_filled]
        )
        # pit node (1, 3) is elevated first above the pass node (3, 1)
        # then sink node (2, 2) is elevated above the pit node
        assert new_elevation[1, 3] > new_elevation[3, 1]
        assert new_elevation[2, 2] > new_elevation[1, 3]

        # only pit node is updated
        assert graph.impl().receivers[8, 0] == 16
        assert graph.impl().receivers_distance[8, 0] > 1e99  # infinity
        assert graph.impl().receivers[12, 0] == 8
        assert graph.impl().receivers_distance[12, 0] == np.sqrt(2)

    @pytest.mark.parametrize("mst_method", [MSTMethod.KRUSKAL, MSTMethod.BORUVKA])
    def test_resolve_carve(self, grid, elevation, mst_method):

        resolver = MSTSinkResolver(mst_method, MSTRouteMethod.CARVE)
        graph = FlowGraph(grid, SingleFlowRouter(), resolver)

        new_elevation = graph.update_routes(elevation)

        # filled elevation (tiny slope)
        mask_filled = np.zeros_like(elevation, dtype=bool)
        mask_filled[[1, 2], [3, 2]] = True
        np.testing.assert_array_equal(
            new_elevation[~mask_filled], elevation[~mask_filled]
        )
        assert new_elevation[2, 2] > new_elevation[3, 1]
        assert new_elevation[1, 3] > new_elevation[2, 2]

        # resolved flow
        assert graph.impl().receivers[8, 0] == 12
        assert graph.impl().receivers_distance[8] == np.sqrt(2)
        assert graph.impl().receivers[12, 0] == 16
        assert graph.impl().receivers_distance[12, 0] == np.sqrt(2)


@pytest.mark.parametrize(
    "resolver",
    [
        PFloodSinkResolver(),
        MSTSinkResolver(MSTMethod.KRUSKAL, MSTRouteMethod.BASIC),
        MSTSinkResolver(MSTMethod.KRUSKAL, MSTRouteMethod.CARVE),
        MSTSinkResolver(MSTMethod.BORUVKA, MSTRouteMethod.BASIC),
        MSTSinkResolver(MSTMethod.BORUVKA, MSTRouteMethod.CARVE),
    ],
)
def test_conservation_drainage_area(resolver):
    # High-level test for sink resolvers that checks the conservation of
    # drainage area (sum of drainage area at domain fixed bounaries = total
    # domain area)

    shape = (101, 101)
    bs = RasterBoundaryStatus(NodeStatus.FIXED_VALUE_BOUNDARY)
    grid = RasterGrid([101, 101], [1, 1], bs, [])

    elevation = np.random.uniform(size=shape)

    graph = FlowGraph(grid, SingleFlowRouter(), resolver)
    graph.update_routes(elevation)
    drainage_area = graph.accumulate(1.0)

    actual = drainage_area[1:-1, [0, -1]].sum() + drainage_area[[0, -1]].sum()
    expected = elevation.size

    assert actual == expected


@pytest.mark.parametrize(
    "resolver",
    [
        PFloodSinkResolver(),
        MSTSinkResolver(MSTMethod.KRUSKAL, MSTRouteMethod.BASIC),
        MSTSinkResolver(MSTMethod.KRUSKAL, MSTRouteMethod.CARVE),
        MSTSinkResolver(MSTMethod.BORUVKA, MSTRouteMethod.BASIC),
        MSTSinkResolver(MSTMethod.BORUVKA, MSTRouteMethod.CARVE),
    ],
)
def test_nb_of_basins(resolver):
    # High-level test for sink resolvers that checks that the total number of
    # basins equals the number of (fixed value) boundary nodes at grid borders
    # (i.e., only outer basins)

    shape = (101, 101)
    bs = RasterBoundaryStatus(NodeStatus.FIXED_VALUE_BOUNDARY)
    grid = RasterGrid([101, 101], [1, 1], bs, [])

    elevation = np.random.uniform(size=shape)

    graph = FlowGraph(grid, SingleFlowRouter(), resolver)
    graph.update_routes(elevation)
    basins = graph.basins()
    nb_basins = basins.max() + 1

    nb_border_nodes = (shape[0] + shape[1]) * 2 - 4

    assert nb_basins == nb_border_nodes
