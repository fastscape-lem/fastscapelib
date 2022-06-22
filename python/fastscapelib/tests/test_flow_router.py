import numpy as np
import numpy.testing as npt
import pytest

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


class TestDummyFlowRouter:
    def test___init__(self):
        dummy_router = DummyFlowRouter()

        with pytest.raises(TypeError):
            DummyFlowRouter(1.0)


class TestSingleFlowRouter:
    @classmethod
    def setup_class(cls):
        profile_grid = ProfileGrid(8, 2.2, [NodeStatus.FIXED_VALUE_BOUNDARY] * 2, [])
        cls.profile_flow_graph = FlowGraph(
            profile_grid, SingleFlowRouter(), NoSinkResolver()
        )
        cls.profile_elevation = np.r_[0.82, 0.16, 0.14, 0.20, 0.71, 0.97, 0.41, 0.09]
        cls.result_profile_elevation = cls.profile_flow_graph.update_routes(
            cls.profile_elevation
        )

        raster_grid = RasterGrid(
            [4, 4],
            [1.1, 1.2],
            RasterBoundaryStatus(NodeStatus.FIXED_VALUE_BOUNDARY),
            [],
        )
        cls.raster_flow_graph = FlowGraph(
            raster_grid, SingleFlowRouter(), NoSinkResolver()
        )
        cls.raster_elevation = np.array(
            [
                [0.82, 0.16, 0.14, 0.20],
                [0.71, 0.97, 0.41, 0.09],
                [0.49, 0.01, 0.19, 0.38],
                [0.29, 0.82, 0.09, 0.88],
            ]
        )
        cls.result_raster_elevation = cls.raster_flow_graph.update_routes(
            cls.raster_elevation
        )

    def test___init__(self):
        single_router = SingleFlowRouter()

        with pytest.raises(TypeError):
            SingleFlowRouter(1.0)

    def test_receivers(self):
        npt.assert_equal(
            self.profile_flow_graph.receivers()[:, 0], np.r_[1, 2, 2, 2, 3, 6, 7, 7]
        )

        npt.assert_equal(
            self.raster_flow_graph.receivers()[:, 0],
            np.r_[1, 2, 7, 7, 9, 9, 7, 7, 9, 9, 9, 7, 9, 9, 9, 14],
        )

    def test_receivers_count(self):
        npt.assert_equal(self.profile_flow_graph.receivers_count(), np.ones(8))

        npt.assert_equal(self.raster_flow_graph.receivers_count(), np.ones(16))

    def test_receivers_distance(self):
        npt.assert_equal(
            self.profile_flow_graph.receivers_distance()[:, 0],
            np.r_[1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0] * 2.2,
        )

        dia = np.sqrt(1.1**2 + 1.2**2)
        npt.assert_equal(
            self.raster_flow_graph.receivers_distance()[:, 0],
            np.r_[
                1.2,
                1.2,
                dia,
                1.1,
                dia,
                1.1,
                1.2,
                0.0,
                1.2,
                0.0,
                1.2,
                1.1,
                dia,
                1.1,
                dia,
                1.2,
            ],
        )

    def test_receivers_weight(self):
        npt.assert_equal(self.profile_flow_graph.receivers_weight()[:, 0], np.ones(8))
        npt.assert_equal(self.profile_flow_graph.receivers_weight()[:, 1], np.zeros(8))

        npt.assert_equal(self.raster_flow_graph.receivers_weight()[:, 0], np.ones(16))
        npt.assert_equal(
            self.raster_flow_graph.receivers_weight()[:, 1:], np.zeros((16, 7))
        )

    def test_donors(self):
        m = np.iinfo(np.uint64).max
        expected_donors = np.full((8, 3), m)
        expected_donors[1, 0] = 0
        expected_donors[2, :] = np.r_[1, 2, 3]
        expected_donors[3, 0] = 4
        expected_donors[6, 0] = 5
        expected_donors[7, :2] = np.r_[6, 7]

        npt.assert_equal(self.profile_flow_graph.donors(), expected_donors)

        expected_donors = np.full((16, 9), m)
        expected_donors[1, 0] = 0
        expected_donors[2, 0] = 1
        expected_donors[7, :5] = np.r_[2, 3, 6, 7, 11]
        expected_donors[9, :8] = np.r_[4, 5, 8, 9, 10, 12, 13, 14]
        expected_donors[14, 0] = 15

        npt.assert_equal(self.raster_flow_graph.donors(), expected_donors)

    def test_donors_count(self):
        npt.assert_equal(
            self.profile_flow_graph.donors_count(), np.r_[0, 1, 3, 1, 0, 0, 1, 2]
        )

        npt.assert_equal(
            self.raster_flow_graph.donors_count(),
            np.r_[0, 1, 1, 0, 0, 0, 0, 5, 0, 8, 0, 0, 0, 0, 1, 0],
        )

    def test_dfs_stack(self):
        npt.assert_equal(
            self.profile_flow_graph.dfs_stack(), np.r_[2, 1, 0, 3, 4, 7, 6, 5]
        )

        npt.assert_equal(
            self.raster_flow_graph.dfs_stack(),
            np.r_[7, 2, 1, 0, 3, 6, 11, 9, 4, 5, 8, 10, 12, 13, 14, 15],
        )


class TestMultipleFlowRouter:
    def test___init__(self):
        multiple_router = MultipleFlowRouter(1.0, 1.5)

        with pytest.raises(TypeError):
            MultipleFlowRouter(1.0)

        with pytest.raises(TypeError):
            MultipleFlowRouter(1.0, "a")
