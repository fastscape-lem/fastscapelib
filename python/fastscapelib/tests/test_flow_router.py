import numpy as np
import numpy.testing as npt
import pytest

from fastscapelib.flow import (
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


class TestSingleFlowRouter:
    @classmethod
    def setup_class(cls):
        profile_grid = ProfileGrid(8, 2.2, [NodeStatus.FIXED_VALUE_BOUNDARY] * 2, [])
        cls.profile_flow_graph = FlowGraph(
            profile_grid, SingleFlowRouter(), NoSinkResolver()
        )
        cls.profile_elevation = np.array([0.0, 0.2, 0.1, 0.2, 0.4, 0.6, 0.3, 0.0])
        cls.result_profile_elevation = cls.profile_flow_graph.update_routes(
            cls.profile_elevation
        )

        # bottom border base level
        bottom_base_level = [
            NodeStatus.CORE,
            NodeStatus.CORE,
            NodeStatus.CORE,
            NodeStatus.FIXED_VALUE_BOUNDARY,
        ]

        raster_grid = RasterGrid(
            [4, 4],
            [1.0, 1.0],
            RasterBoundaryStatus(bottom_base_level),
            [],
        )
        cls.raster_flow_graph = FlowGraph(
            raster_grid, SingleFlowRouter(), NoSinkResolver()
        )

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

    def test___init__(self):
        SingleFlowRouter()

        # no parameter
        with pytest.raises(TypeError):
            SingleFlowRouter(1.0)

    def test_receivers(self):
        npt.assert_equal(
            self.profile_flow_graph.impl().receivers()[:, 0],
            np.array([0, 0, 2, 2, 3, 6, 7, 7]),
        )

        npt.assert_equal(
            self.raster_flow_graph.impl().receivers()[:, 0],
            np.array([4, 5, 6, 7, 8, 9, 10, 11, 13, 13, 13, 15, 12, 13, 14, 15]),
        )

    def test_receivers_count(self):
        npt.assert_equal(self.profile_flow_graph.impl().receivers_count(), np.ones(8))

        npt.assert_equal(self.raster_flow_graph.impl().receivers_count(), np.ones(16))

    def test_receivers_distance(self):
        npt.assert_equal(
            self.profile_flow_graph.impl().receivers_distance()[:, 0],
            np.array([0, 1, 0, 1, 1, 1, 1, 0]) * 2.2,
        )

        dia = np.sqrt(2)
        npt.assert_equal(
            self.raster_flow_graph.impl().receivers_distance()[:, 0],
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

    def test_receivers_weight(self):
        npt.assert_equal(
            self.profile_flow_graph.impl().receivers_weight()[:, 0], np.ones(8)
        )
        npt.assert_equal(
            self.raster_flow_graph.impl().receivers_weight()[:, 0], np.ones(16)
        )

    def test_donors(self):
        m = np.iinfo(np.uint64).max
        expected_donors = (
            np.array(
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
            ),
        )
        npt.assert_equal(
            self.profile_flow_graph.impl().donors(), np.squeeze(expected_donors)
        )

        expected_donors = (
            np.array(
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
            ),
        )
        npt.assert_equal(
            self.raster_flow_graph.impl().donors(), np.squeeze(expected_donors)
        )

    def test_donors_count(self):
        npt.assert_equal(
            self.profile_flow_graph.impl().donors_count(),
            np.array([1, 0, 2, 1, 0, 0, 1, 1]),
        )

        npt.assert_equal(
            self.raster_flow_graph.impl().donors_count(),
            np.array([0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 3, 0, 1]),
        )

    def test_dfs_indices(self):
        npt.assert_equal(
            self.profile_flow_graph.impl().dfs_indices(),
            np.array([0, 1, 2, 3, 4, 7, 6, 5]),
        )

        npt.assert_equal(
            self.raster_flow_graph.impl().dfs_indices(),
            np.array([12, 13, 8, 9, 10, 6, 2, 5, 1, 4, 0, 14, 15, 11, 7, 3]),
        )


class TestMultipleFlowRouter:
    def test___init__(self):
        MultipleFlowRouter(1.0, 1.5)

        with pytest.raises(TypeError):
            MultipleFlowRouter(1.0)

        with pytest.raises(TypeError):
            MultipleFlowRouter(1.0, "a")
