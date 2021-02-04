import pytest
import numpy as np
import numpy.testing as npt

from fastscapelib.grid import NodeStatus, RasterBoundaryStatus, RasterNode, RasterNeighbor, RasterGrid


class TestRasterBoundaryStatus():

    def setup_method(self, method):
        self.bs1 = RasterBoundaryStatus([NodeStatus.CORE, NodeStatus.FIXED_VALUE_BOUNDARY, NodeStatus.CORE, NodeStatus.FIXED_VALUE_BOUNDARY])
        self.bs2 = RasterBoundaryStatus(NodeStatus.CORE)
        self.bs3 = RasterBoundaryStatus(NodeStatus.LOOPED_BOUNDARY)

    def test__init__(self):
        # call the setup_method
        pass

    def test_left(self):
        assert self.bs1.left == NodeStatus.CORE
        assert self.bs2.left == NodeStatus.CORE
        assert self.bs3.left == NodeStatus.LOOPED_BOUNDARY

    def test_right(self):
        assert self.bs1.right == NodeStatus.FIXED_VALUE_BOUNDARY
        assert self.bs2.right == NodeStatus.CORE
        assert self.bs3.right == NodeStatus.LOOPED_BOUNDARY

    def test_bottom(self):
        assert self.bs1.bottom == NodeStatus.FIXED_VALUE_BOUNDARY
        assert self.bs2.bottom == NodeStatus.CORE
        assert self.bs3.bottom == NodeStatus.LOOPED_BOUNDARY

    def test_top(self):
        assert self.bs1.top == NodeStatus.CORE
        assert self.bs2.top == NodeStatus.CORE
        assert self.bs3.top == NodeStatus.LOOPED_BOUNDARY

    def test_is_horizontal_looped(self):
        assert not self.bs1.is_horizontal_looped
        assert not self.bs2.is_horizontal_looped
        assert self.bs3.is_horizontal_looped

    def test_is_vertical_looped(self):
        assert not self.bs1.is_vertical_looped
        assert not self.bs2.is_vertical_looped
        assert self.bs3.is_vertical_looped


class TestRasterNeighbor():

    def setup_method(self, method):
        self.n = RasterNeighbor(5, 0, 5, 1.35, NodeStatus.CORE)

    def test__init__(self):
        # call the setup_method
        assert self.n == RasterNeighbor(5, 0, 5, 1.35, NodeStatus.CORE)
        pass

    def test_flatten_idx(self):
        assert self.n.flatten_idx == 5

        self.n.flatten_idx = 8
        assert self.n.flatten_idx == 8

    def test_row(self):
        assert self.n.row == 0

        self.n.row = 3
        assert self.n.row == 3

    def test_col(self):
        assert self.n.col == 5

        self.n.col = 8
        assert self.n.col == 8

    def test_distance(self):
        assert self.n.distance == 1.35

        self.n.distance = 1.38
        assert self.n.distance == 1.38

    def test_status(self):
        assert self.n.status == NodeStatus.CORE

        self.n.status = NodeStatus.FIXED_VALUE_BOUNDARY
        assert self.n.status == NodeStatus.FIXED_VALUE_BOUNDARY


class TestRasterGrid():

    def setup_method(self, method):
        self.bs = bs = RasterBoundaryStatus(NodeStatus.FIXED_VALUE_BOUNDARY)
        self.g = RasterGrid([5, 10], [2.2, 2.4], bs, [RasterNode(0, 5, NodeStatus.FIXED_VALUE_BOUNDARY)])

    def test___init__(self):
        g1 = RasterGrid([10, 10], [2.3, 2.1], self.bs, [RasterNode(0, 5, NodeStatus.FIXED_VALUE_BOUNDARY)])
        assert g1.size == 100
        npt.assert_almost_equal(g1.spacing, [2.3, 2.1])
        npt.assert_almost_equal(g1.length, [20.7, 18.9])

        with pytest.raises(IndexError):
            RasterGrid([5, 10], [2.2, 2.4], self.bs, [RasterNode(20, 255, NodeStatus.FIXED_VALUE_BOUNDARY)])

    def test_from_length(self):
        g = RasterGrid.from_length([11, 11], np.r_[23, 21], self.bs, [RasterNode(0, 5, NodeStatus.FIXED_VALUE_BOUNDARY)])
        assert g.size == 121
        npt.assert_almost_equal(g.spacing, [2.3, 2.1])
        npt.assert_almost_equal(g.length, [23., 21.])

    def test_size(self):
        assert self.g.size == 50

    def test_shape(self):
        npt.assert_almost_equal(self.g.shape, [5, 10])

    def test_spacing(self):
        npt.assert_almost_equal(self.g.spacing, [2.2, 2.4])
