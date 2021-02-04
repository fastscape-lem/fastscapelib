import pytest
import numpy as np

from fastscapelib.algo import erode_stream_power_d, erode_stream_power_var_d
from fastscapelib.grid import ProfileGrid, RasterGrid, Node, NodeStatus, RasterBoundaryStatus, RasterNode
from fastscapelib.flow_graph import FlowGraph, DummyFlowRouter, MultipleFlowRouter, SingleFlowRouter


class TestErodeStreamPower():

    @pytest.mark.parametrize("func,k", [(erode_stream_power_d, 1e-3),
                                        (erode_stream_power_var_d, np.full(4, 1e-3))])
    def test_profile_grid(self, func, k):
        spacing = 300.
        grid = ProfileGrid(4, 300, [NodeStatus.FIXED_VALUE_BOUNDARY]*2, [])
        flow_graph = FlowGraph(grid, SingleFlowRouter())
        
        h = 1.
        elevation = np.array([0., h, h, 0.], dtype='d')
        graph_elevation = flow_graph.update_routes(elevation)

        drainage_area = flow_graph.accumulate(1.)
        erosion = np.zeros_like(elevation)
        m_exp = 0.5
        n_exp = 1.

        dt = 1.   # use small time step (compare with explicit scheme)
        tolerance = 1e-3

        n_corr = func(erosion, elevation, drainage_area, flow_graph,
                      k, m_exp, n_exp, dt, tolerance)

        slope = h / spacing
        a = spacing
        k_coef = 1e-3
        err = dt * k_coef * a**m_exp * slope**n_exp
        expected_erosion = np.array([0., err, err, 0.], dtype='d')

        np.testing.assert_allclose(erosion, expected_erosion, atol=1e-5)
        assert n_corr == 0

    @pytest.mark.parametrize("func,k", [(erode_stream_power_d, 1e-3),
                                        (erode_stream_power_var_d, np.full((2, 2), 1e-3))])
    def test_raster_grid(self, func, k):
        # Test on a tiny (2x2) 2-d square grid with a planar surface
        # tilted in y (rows) and with all outlets on the 1st row.
        spacing = 300.
        grid = RasterGrid([2, 2], [spacing, spacing], RasterBoundaryStatus(NodeStatus.FIXED_VALUE_BOUNDARY), [])
        flow_graph = FlowGraph(grid, SingleFlowRouter())

        h = 1.
        elevation = np.array([[0., 0.], [h, h]], dtype='d')
        graph_elevation = flow_graph.update_routes(elevation)

        drainage_area = flow_graph.accumulate(1.)
        erosion = np.zeros_like(elevation)
        m_exp = 0.5
        n_exp = 1.

        dt = 1   # use small time step (compare with explicit scheme)
        tolerance = 1e-3

        n_corr = func(erosion, elevation, drainage_area, flow_graph,
                      k, m_exp, n_exp, dt, tolerance)

        slope = h / spacing
        a = spacing**2
        k_coef = 1e-3
        err = dt * k_coef * a**m_exp * slope**n_exp
        expected_erosion = np.array([[0., 0.], [err, err]], dtype='d')

        np.testing.assert_allclose(erosion, expected_erosion, atol=1e-5)
        assert n_corr == 0