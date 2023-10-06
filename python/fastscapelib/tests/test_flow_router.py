import numpy as np
import numpy.testing as npt
import pytest

from fastscapelib.flow import (
    FlowDirection,
    FlowGraph,
    FlowOperator,
    MultiFlowRouter,
    PFloodSinkResolver,
    SingleFlowRouter,
)
from fastscapelib.grid import NodeStatus, ProfileGrid, RasterGrid


class TestSingleFlowRouter:
    profile_flow_graph: FlowGraph
    profile_elevation: np.ndarray
    result_profile_elevation: np.ndarray
    raster_flow_graph: FlowGraph
    raster_elevation: np.ndarray
    result_raster_elevation: np.ndarray

    @classmethod
    def setup_class(cls) -> None:
        profile_grid = ProfileGrid(8, 2.2, NodeStatus.FIXED_VALUE)
        cls.profile_flow_graph = FlowGraph(profile_grid, [SingleFlowRouter()])
        cls.profile_elevation = np.array([0.0, 0.2, 0.1, 0.2, 0.4, 0.6, 0.3, 0.0])
        cls.result_profile_elevation = cls.profile_flow_graph.update_routes(
            cls.profile_elevation
        )

        # bottom border base level
        bottom_base_level = [
            NodeStatus.CORE,
            NodeStatus.CORE,
            NodeStatus.CORE,
            NodeStatus.FIXED_VALUE,
        ]

        raster_grid = RasterGrid([4, 4], [1.0, 1.0], bottom_base_level)

        cls.raster_flow_graph = FlowGraph(raster_grid, [SingleFlowRouter()])

        # planar surface tilted along the y-axis + small carved channel
        cls.raster_elevation = np.array(
            [
                [0.6, 0.6, 0.6, 0.6],
                [0.4, 0.4, 0.4, 0.4],
                [0.2, 0.2, 0.2, 0.2],
                [0.1, 0.0, 0.1, 0.1],
            ]
        )
        cls.result_raster_elevation = cls.raster_flow_graph.update_routes(
            cls.raster_elevation
        )

    def test_constructor(self) -> None:
        router = SingleFlowRouter()
        assert isinstance(router, FlowOperator)
        assert router.name == "single_flow_router"
        assert repr(router) == "SingleFlowRouter"

    def test_class_attrs(self) -> None:
        assert SingleFlowRouter.graph_updated is True
        assert SingleFlowRouter.elevation_updated is False
        assert SingleFlowRouter.in_flowdir == FlowDirection.UNDEFINED
        assert SingleFlowRouter.out_flowdir == FlowDirection.SINGLE

    def test_receivers(self) -> None:
        npt.assert_equal(
            self.profile_flow_graph.impl().receivers[:, 0],
            np.array([0, 0, 2, 2, 3, 6, 7, 7]),
        )

        npt.assert_equal(
            self.raster_flow_graph.impl().receivers[:, 0],
            np.array([4, 5, 6, 7, 8, 9, 10, 11, 13, 13, 13, 15, 12, 13, 14, 15]),
        )

    def test_receivers_count(self) -> None:
        npt.assert_equal(self.profile_flow_graph.impl().receivers_count, np.ones(8))

        npt.assert_equal(self.raster_flow_graph.impl().receivers_count, np.ones(16))

    def test_receivers_distance(self) -> None:
        npt.assert_equal(
            self.profile_flow_graph.impl().receivers_distance[:, 0],
            np.array([0, 1, 0, 1, 1, 1, 1, 0]) * 2.2,
        )

        dia = np.sqrt(2)
        npt.assert_equal(
            self.raster_flow_graph.impl().receivers_distance[:, 0],
            np.array(
                [
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    dia,
                    1.0,
                    dia,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ]
            ),
        )

    def test_receivers_weight(self) -> None:
        npt.assert_equal(
            self.profile_flow_graph.impl().receivers_weight[:, 0], np.ones(8)
        )
        npt.assert_equal(
            self.raster_flow_graph.impl().receivers_weight[:, 0], np.ones(16)
        )

    def test_donors(self) -> None:
        actual = self.profile_flow_graph.impl().donors
        m = np.iinfo(actual.dtype).max
        expected = np.array(
            [
                [1, m, m],
                [m, m, m],
                [2, 3, m],
                [4, m, m],
                [m, m, m],
                [m, m, m],
                [5, m, m],
                [6, m, m],
            ]
        )
        npt.assert_equal(actual, expected)

        actual = self.raster_flow_graph.impl().donors
        m = np.iinfo(actual.dtype).max
        expected = np.array(
            [
                [m, m, m, m, m, m, m, m, m],
                [m, m, m, m, m, m, m, m, m],
                [m, m, m, m, m, m, m, m, m],
                [m, m, m, m, m, m, m, m, m],
                [0, m, m, m, m, m, m, m, m],
                [1, m, m, m, m, m, m, m, m],
                [2, m, m, m, m, m, m, m, m],
                [3, m, m, m, m, m, m, m, m],
                [4, m, m, m, m, m, m, m, m],
                [5, m, m, m, m, m, m, m, m],
                [6, m, m, m, m, m, m, m, m],
                [7, m, m, m, m, m, m, m, m],
                [m, m, m, m, m, m, m, m, m],
                [8, 9, 10, m, m, m, m, m, m],
                [m, m, m, m, m, m, m, m, m],
                [11, m, m, m, m, m, m, m, m],
            ]
        )
        npt.assert_equal(actual, expected)

    def test_donors_count(self) -> None:
        npt.assert_equal(
            self.profile_flow_graph.impl().donors_count,
            np.array([1, 0, 2, 1, 0, 0, 1, 1]),
        )

        npt.assert_equal(
            self.raster_flow_graph.impl().donors_count,
            np.array([0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 3, 0, 1]),
        )

    def test_dfs_indices(self) -> None:
        npt.assert_equal(
            self.profile_flow_graph.impl().dfs_indices,
            np.array([0, 1, 2, 3, 4, 7, 6, 5]),
        )

        npt.assert_equal(
            self.raster_flow_graph.impl().dfs_indices,
            np.array([12, 13, 8, 9, 10, 6, 2, 5, 1, 4, 0, 14, 15, 11, 7, 3]),
        )


class TestMultiFlowRouter:
    def test_constructor(self) -> None:
        router = MultiFlowRouter()
        assert isinstance(router, FlowOperator)
        assert router.name == "multi_flow_router"
        assert repr(router).startswith("MultiFlowRouter (slope_exp=1.0")
        assert router.slope_exp == 1.0

        router = MultiFlowRouter(0.0)
        assert router.slope_exp == 0.0

        router.slope_exp = 0.5
        assert router.slope_exp == 0.5

    def test_class_attrs(self) -> None:
        assert MultiFlowRouter.graph_updated is True
        assert MultiFlowRouter.elevation_updated is False
        assert MultiFlowRouter.in_flowdir == FlowDirection.UNDEFINED
        assert MultiFlowRouter.out_flowdir == FlowDirection.MULTI

    def test_graph_topology(self) -> None:
        # test on a 3x3 tiny grid with base levels at all border nodes
        grid = RasterGrid([3, 3], [1.0, 1.0], NodeStatus.FIXED_VALUE)

        diag = np.sqrt(2)
        elevation = np.array([[0.0, 1.0, 0.0], [1.0, 2.0, 1.0], [0.0, 1.0, 0.0]])
        router = MultiFlowRouter(0.0)
        flow_graph = FlowGraph(grid, [router])
        flow_graph.update_routes(elevation)

        npt.assert_array_equal(
            flow_graph.impl().receivers_count,
            [1, 1, 1, 1, 8, 1, 1, 1, 1],
        )
        npt.assert_array_equal(
            flow_graph.impl().receivers[4],
            [0, 1, 2, 3, 5, 6, 7, 8],
        )
        npt.assert_array_almost_equal(
            flow_graph.impl().receivers_distance[4],
            [diag, 1.0, diag, 1.0, 1.0, diag, 1.0, diag],
        )
        npt.assert_array_equal(
            flow_graph.impl().donors_count,
            [1, 1, 1, 1, 0, 1, 1, 1, 1],
        )
        npt.assert_array_equal(
            flow_graph.impl().donors[[0, 1, 2, 3, 5, 6, 7, 8], 0].ravel(),
            [4, 4, 4, 4, 4, 4, 4, 4],
        )

    def test_mask(self) -> None:
        # see test_graph_topology for test setup
        grid = RasterGrid([3, 3], [1.0, 1.0], NodeStatus.FIXED_VALUE)

        elevation = np.array([[0.0, 1.0, 0.0], [1.0, 2.0, 1.0], [0.0, 1.0, 0.0]])

        router = MultiFlowRouter(0.0)
        flow_graph = FlowGraph(grid, [router])

        # masked center node -> no graph edge (each node has itself as receiver)
        mask = np.array(
            [[False, False, False], [False, True, False], [False, False, False]]
        )
        flow_graph.mask = mask
        flow_graph.update_routes(elevation)
        npt.assert_array_equal(flow_graph.impl().receivers_count, np.ones(9, dtype=int))

        # one masked edge node -> this masked node not in center node receivers
        mask = np.array(
            [[False, False, False], [True, False, False], [False, False, False]]
        )
        flow_graph.mask = mask
        flow_graph.update_routes(elevation)
        npt.assert_array_equal(
            flow_graph.impl().receivers[4][:-1],
            [0, 1, 2, 5, 6, 7, 8],
        )

    @pytest.mark.parametrize(
        "slope_exp, weights",
        [
            (0.0, [1 / 8, 1 / 8, 1 / 8, 1 / 8, 1 / 8, 1 / 8, 1 / 8, 1 / 8]),
            (2.0, [2 / 12, 1 / 12, 2 / 12, 1 / 12, 1 / 12, 2 / 12, 1 / 12, 2 / 12]),
        ],
    )
    def test_flow_partitions(self, slope_exp, weights) -> None:
        # test on a 3x3 tiny grid with base levels at all border nodes
        grid = RasterGrid([3, 3], [1.0, 1.0], NodeStatus.FIXED_VALUE)

        sqrt2 = np.sqrt(2)
        elevation = np.array(
            [[0.0, sqrt2, 0.0], [sqrt2, sqrt2 * 2, sqrt2], [0.0, sqrt2, 0.0]]
        )
        router = MultiFlowRouter(slope_exp)
        flow_graph = FlowGraph(grid, [router])
        flow_graph.update_routes(elevation)

        npt.assert_allclose(
            flow_graph.impl().receivers_weight[4],
            weights,
        )


@pytest.mark.parametrize(
    "router", [SingleFlowRouter(), MultiFlowRouter(0.0), MultiFlowRouter(2.0)]
)
def test_conservation_area(router) -> None:
    # High level test: conservative flow routing
    nrows = 10
    ncols = 8

    # bottom border base level
    bottom_base_level = [
        NodeStatus.CORE,
        NodeStatus.CORE,
        NodeStatus.CORE,
        NodeStatus.FIXED_VALUE,
    ]
    grid = RasterGrid([nrows, ncols], [1.0, 1.0], bottom_base_level)

    # avoid closed depressions (all flow must reach bottom border nodes)
    flow_graph = FlowGraph(grid, [PFloodSinkResolver(), router])

    # planar surface tilted along the y-axis + random perturbations
    elevation = np.random.uniform(size=grid.shape) + np.arange(nrows)[:, None] * 2

    flow_graph.update_routes(elevation)
    drainage_area = flow_graph.accumulate(1.0)

    # assumes grid cell uniform area is 1
    assert abs(np.sum(drainage_area[-1]) - grid.size) < 1e-5


@pytest.mark.parametrize(
    "router", [SingleFlowRouter(), MultiFlowRouter(0.0), MultiFlowRouter(2.0)]
)
def test_monotonic_dfs(router) -> None:
    # High level test: monotonic elevation for dfs indices
    nrows = 10
    ncols = 8

    grid = RasterGrid([nrows, ncols], [1.0, 1.0], NodeStatus.FIXED_VALUE)

    # this test requires no closed depressions (fill them with tiny slope)
    flow_graph = FlowGraph(grid, [PFloodSinkResolver(), router])

    elevation = np.random.uniform(size=grid.shape)
    filled_elevation = flow_graph.update_routes(elevation)
    filled_elevation_flat = filled_elevation.ravel()

    receivers_count = flow_graph.impl().receivers_count
    receivers = flow_graph.impl().receivers
    dfs_indices = flow_graph.impl().dfs_indices

    # traverse the graph, check that elevation increases if
    # previous node is a direct receiver of current visited node
    # (dfs_indices always in bottom->up direction)
    for i, inode in enumerate(dfs_indices[1:]):
        iprev = dfs_indices[i - 1]
        if iprev in receivers[inode][: receivers_count[inode]]:
            assert filled_elevation_flat[inode] >= filled_elevation_flat[iprev]
