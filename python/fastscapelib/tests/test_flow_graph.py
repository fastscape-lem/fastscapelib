import pytest
import numpy as np
import numpy.testing as npt

from fastscapelib.grid import ProfileGrid, RasterGrid, Node, NodeStatus, RasterBoundaryStatus, RasterNode
from fastscapelib.flow_graph import FlowGraph, DummyFlowRouter, MultipleFlowRouter, SingleFlowRouter


class TestFlowGraph():

    def test___init__(self): 
        profile_grid = ProfileGrid(8, 2.2, [NodeStatus.FIXED_VALUE_BOUNDARY]*2, [])
        raster_grid = RasterGrid([5, 10], [2.2, 2.4], RasterBoundaryStatus(NodeStatus.FIXED_VALUE_BOUNDARY), [])

        FlowGraph(profile_grid, DummyFlowRouter())
        FlowGraph(profile_grid, MultipleFlowRouter(1., 1.1))
        FlowGraph(profile_grid, SingleFlowRouter())

    def test_update_routes(self):
        grid = ProfileGrid(8, 2.2, [NodeStatus.FIXED_VALUE_BOUNDARY]*2, [])
        flow_graph = FlowGraph(grid, SingleFlowRouter())
        elevation = np.r_[0.82,  0.16,  0.14,  0.20, 0.71,  0.97,  0.41,  0.09]

        graph_elevation = flow_graph.update_routes(elevation)

        npt.assert_equal(flow_graph.receivers()[:, 0], np.r_[1, 2, 2, 2, 3, 6, 7, 7])
        npt.assert_equal(flow_graph.receivers_count(), np.ones(elevation.size))
        npt.assert_equal(flow_graph.receivers_weight()[:, 0], np.ones(elevation.size))
        npt.assert_equal(flow_graph.receivers_weight()[:, 1], np.zeros(elevation.size))

        m = np.iinfo(np.uint64).max
        npt.assert_equal(flow_graph.donors(), np.array([[m, m, m],
                                                        [0, m, m],
                                                        [1, 2, 3],
                                                        [4, m, m],
                                                        [m, m, m],
                                                        [m, m, m],
                                                        [5, m, m],
                                                        [6, 7, m]]))
        npt.assert_equal(flow_graph.donors_count(), np.r_[0, 1, 3, 1, 0, 0, 1, 2])

        npt.assert_equal(graph_elevation, elevation)

    def test_accumulate(self):
        grid = ProfileGrid(8, 2.2, [NodeStatus.FIXED_VALUE_BOUNDARY]*2, [])
        flow_graph = FlowGraph(grid, SingleFlowRouter())
        elevation = np.r_[0.82,  0.16,  0.14,  0.20, 0.71,  0.97,  0.41,  0.09]

        graph_elevation = flow_graph.update_routes(elevation)

        npt.assert_almost_equal(flow_graph.accumulate(np.ones(elevation.shape)),
                                np.r_[2.2, 4.4, 11., 4.4, 2.2, 2.2, 4.4, 6.6])

        grid = RasterGrid([4, 4], [1.1, 1.2], RasterBoundaryStatus(NodeStatus.FIXED_VALUE_BOUNDARY), [])
        flow_graph = FlowGraph(grid, SingleFlowRouter())
        elevation = np.array([[0.82,  0.16,  0.14,  0.20],
                              [0.71,  0.97,  0.41,  0.09],
                              [0.49,  0.01,  0.19,  0.38],
                              [0.29,  0.82,  0.09,  0.88]])

        graph_elevation = flow_graph.update_routes(elevation)

        expected = np.array([[1.32,  2.64, 3.96, 1.32],
                             [1.32,  1.32, 1.32, 9.24],
                             [1.32, 11.88, 1.32, 1.32],
                             [1.32,  1.32, 2.64, 1.32]])
        npt.assert_almost_equal(flow_graph.accumulate(np.ones(elevation.shape)),
                                expected)

        npt.assert_almost_equal(flow_graph.accumulate(1.), expected)

        npt.assert_almost_equal(flow_graph.accumulate(5.), 5 * expected)

        data = np.array([[1.1, 1.0, 1.1, 1.0],
                         [1.1, 1.0, 1.1, 1.0],
                         [1.1, 1.0, 1.1, 1.0],
                         [1.1, 1.0, 1.1, 1.0]])
        expected = np.array([[1.452, 2.772, 4.224,  1.32],
                             [1.452,  1.32, 1.452, 9.636],
                             [1.452, 12.54, 1.452,  1.32],
                             [1.452,  1.32, 2.772,  1.32]])
        npt.assert_almost_equal(flow_graph.accumulate(data),
                                expected)
