import pytest
import numpy as np
import numpy.testing as npt

import fastscapelib
from _fastscapelib_py.grid import NodeStatus, ProfileBoundaryStatus, ProfileGrid, Node, Neighbor


class TestProfileBoundaryStatus():

    def setup_method(self, method):
        self.bs1 = ProfileBoundaryStatus([NodeStatus.CORE, NodeStatus.FIXED_VALUE_BOUNDARY])
        self.bs2 = ProfileBoundaryStatus(NodeStatus.FIXED_VALUE_BOUNDARY, NodeStatus.CORE)
        self.bs3 = ProfileBoundaryStatus(NodeStatus.CORE)
        self.bs4 = ProfileBoundaryStatus(NodeStatus.LOOPED_BOUNDARY)

    def test__init__(self):
        # call the setup_method
        pass

    def test_left(self):
        assert self.bs1.left == NodeStatus.CORE
        assert self.bs2.left == NodeStatus.FIXED_VALUE_BOUNDARY
        assert self.bs3.left == NodeStatus.CORE
        assert self.bs4.left == NodeStatus.LOOPED_BOUNDARY

    def test_right(self):
        assert self.bs1.right == NodeStatus.FIXED_VALUE_BOUNDARY
        assert self.bs2.right == NodeStatus.CORE
        assert self.bs3.right == NodeStatus.CORE
        assert self.bs4.right == NodeStatus.LOOPED_BOUNDARY

    def test_is_horizontal_looped(self):
        assert not self.bs1.is_horizontal_looped
        assert not self.bs2.is_horizontal_looped
        assert not self.bs3.is_horizontal_looped
        assert self.bs4.is_horizontal_looped


class TestNeighbor():

    def setup_method(self, method):
        self.n = Neighbor(5, 1.35, NodeStatus.CORE)

    def test__init__(self):
        # call the setup_method
        assert self.n == Neighbor(5, 1.35, NodeStatus.CORE)
        pass

    def test_idx(self):
        assert self.n.idx == 5

        self.n.idx = 8
        assert self.n.idx == 8

    def test_distance(self):
        assert self.n.distance == 1.35

        self.n.distance = 1.38
        assert self.n.distance == 1.38

    def test_status(self):
        assert self.n.status == NodeStatus.CORE

        self.n.status = NodeStatus.FIXED_VALUE_BOUNDARY
        assert self.n.status == NodeStatus.FIXED_VALUE_BOUNDARY


class TestProfileGrid():

    def setup_method(self, method):
        bs = [NodeStatus.FIXED_VALUE_BOUNDARY]*2
        self.g = ProfileGrid(10, 2.2, bs, [(5, NodeStatus.FIXED_VALUE_BOUNDARY)])

    def test___init__(self):
        bs = [NodeStatus.FIXED_VALUE_BOUNDARY]*2
        g = ProfileGrid(10, 2, bs, [(5, NodeStatus.FIXED_VALUE_BOUNDARY)])
        assert g.size == 10
        assert g.spacing == 2.

        g = ProfileGrid(15, 3., ProfileBoundaryStatus(bs), [Node(5, NodeStatus.FIXED_VALUE_BOUNDARY)])
        assert g.size == 15
        assert g.spacing == 3.

        with pytest.raises(IndexError):
            ProfileGrid(10, 2, bs, [(15, NodeStatus.FIXED_VALUE_BOUNDARY)])

    def test_size(self):
        assert self.g.size == 10

    def test_spacing(self):
        assert self.g.spacing == 2.2

    def test_status_at_nodes(self):
        npt.assert_equal(self.g.status_at_nodes, np.array([1, 0, 0, 0, 0, 1, 0, 0, 0, 1]))

#    def test_neighbors(self):
#        assert self.g.neighbors(0) == [Neighbor(1, 2.2, NodeStatus.CORE)]
#        assert self.g.neighbors(1) == [Neighbor(0, 2.2, NodeStatus.FIXED_VALUE_BOUNDARY), Neighbor(2, 2.2, NodeStatus.CORE)]
#        assert self.g.neighbors(6) == [Neighbor(5, 2.2, NodeStatus.FIXED_VALUE_BOUNDARY), Neighbor(7, 2.2, NodeStatus.CORE)]
#        
#        with pytest.raises(IndexError):
#            self.g.neighbors(11)
