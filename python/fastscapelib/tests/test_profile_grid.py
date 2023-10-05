import numpy as np
import numpy.testing as npt
import pytest

from fastscapelib.grid import (
    Neighbor,
    Node,
    NodeStatus,
    ProfileBoundaryStatus,
    ProfileGrid,
)


class TestProfileBoundaryStatus:
    def setup_method(self) -> None:
        self.bs1 = ProfileBoundaryStatus([NodeStatus.CORE, NodeStatus.FIXED_VALUE])
        self.bs2 = ProfileBoundaryStatus(NodeStatus.FIXED_VALUE, NodeStatus.CORE)
        self.bs3 = ProfileBoundaryStatus(NodeStatus.CORE)
        self.bs4 = ProfileBoundaryStatus(NodeStatus.LOOPED)

    def test__init__(self) -> None:
        # call the setup_method
        pass

    def test_left(self) -> None:
        assert self.bs1.left == NodeStatus.CORE
        assert self.bs2.left == NodeStatus.FIXED_VALUE
        assert self.bs3.left == NodeStatus.CORE
        assert self.bs4.left == NodeStatus.LOOPED

    def test_right(self) -> None:
        assert self.bs1.right == NodeStatus.FIXED_VALUE
        assert self.bs2.right == NodeStatus.CORE
        assert self.bs3.right == NodeStatus.CORE
        assert self.bs4.right == NodeStatus.LOOPED

    def test_is_horizontal_looped(self) -> None:
        assert not self.bs1.is_horizontal_looped
        assert not self.bs2.is_horizontal_looped
        assert not self.bs3.is_horizontal_looped
        assert self.bs4.is_horizontal_looped


class TestNeighbor:
    def setup_method(self) -> None:
        self.n = Neighbor(5, 1.35, NodeStatus.CORE)

    def test__init__(self) -> None:
        # call the setup_method
        assert self.n == Neighbor(5, 1.35, NodeStatus.CORE)
        pass

    def test_idx(self) -> None:
        assert self.n.idx == 5

        self.n.idx = 8
        assert self.n.idx == 8

    def test_distance(self) -> None:
        assert self.n.distance == 1.35

        self.n.distance = 1.38
        assert self.n.distance == 1.38

    def test_status(self) -> None:
        assert self.n.status == NodeStatus.CORE

        self.n.status = NodeStatus.FIXED_VALUE
        assert self.n.status == NodeStatus.FIXED_VALUE


class TestProfileGrid:
    def setup_method(self) -> None:
        self.bs = [NodeStatus.FIXED_VALUE] * 2
        self.g = ProfileGrid(10, 2.2, self.bs, nodes_status={5: NodeStatus.FIXED_VALUE})

    def test_static_props(self) -> None:
        assert ProfileGrid.is_structured is True
        assert ProfileGrid.is_uniform is True
        assert ProfileGrid.n_neighbors_max == 2

    def test_constructor(self) -> None:
        g = ProfileGrid(10, 2, self.bs, nodes_status={5: NodeStatus.FIXED_VALUE})
        assert g.size == 10
        assert g.spacing == 2.0
        assert g.length == 18.0

        g = ProfileGrid(
            15,
            3.0,
            ProfileBoundaryStatus(self.bs),
            {5: NodeStatus.FIXED_VALUE},
        )
        assert g.size == 15
        assert g.spacing == 3.0
        assert g.length == 42.0

        with pytest.raises(IndexError):
            ProfileGrid(10, 2, self.bs, {15: NodeStatus.FIXED_VALUE})

    @pytest.mark.parametrize(
        "bounds_status",
        [
            NodeStatus.FIXED_VALUE,
            [NodeStatus.FIXED_VALUE] * 2,
            ProfileBoundaryStatus(NodeStatus.FIXED_VALUE),
        ],
    )
    def test_constructor_bounds_status(self, bounds_status) -> None:
        g = ProfileGrid(10, 2, bounds_status)
        assert g.nodes_status(0) == g.nodes_status(9) == NodeStatus.FIXED_VALUE

    def test_from_length(self) -> None:
        g = ProfileGrid.from_length(
            11,
            20.0,
            ProfileBoundaryStatus(self.bs),
            {5: NodeStatus.FIXED_VALUE},
        )
        assert g.size == 11
        assert g.shape == [11]
        assert g.spacing == 2.0
        assert g.length == 20.0

    @pytest.mark.parametrize(
        "bounds_status",
        [
            NodeStatus.FIXED_VALUE,
            [NodeStatus.FIXED_VALUE] * 2,
            ProfileBoundaryStatus(NodeStatus.FIXED_VALUE),
        ],
    )
    def test_from_length_bounds_status(self, bounds_status) -> None:
        g = ProfileGrid.from_length(10, 20.0, bounds_status)
        assert g.nodes_status(0) == g.nodes_status(9) == NodeStatus.FIXED_VALUE

    def test_size(self) -> None:
        assert self.g.size == 10

    def test_spacing(self) -> None:
        assert self.g.spacing == 2.2

    def test_shape(self) -> None:
        npt.assert_almost_equal(self.g.shape, np.r_[10])

    def test_length(self) -> None:
        assert self.g.length == 19.8

    def test_nodes_indices(self) -> None:
        npt.assert_equal(self.g.nodes_indices(), np.arange(self.g.size))
        npt.assert_equal(self.g.nodes_indices(NodeStatus.FIXED_VALUE), [0, 5, 9])
        npt.assert_equal(self.g.nodes_indices(NodeStatus.CORE), [1, 2, 3, 4, 6, 7, 8])
        assert not len(self.g.nodes_indices(NodeStatus.FIXED_GRADIENT))

    def test_nodes_status(self) -> None:
        npt.assert_equal(
            self.g.nodes_status(), np.array([1, 0, 0, 0, 0, 1, 0, 0, 0, 1])
        )
        # should return a copy
        assert not np.shares_memory(self.g.nodes_status(), self.g.nodes_status())

        assert self.g.nodes_status(0) == NodeStatus.FIXED_VALUE
        assert self.g.nodes_status(1) == NodeStatus.CORE

    def test_nodes_areas(self) -> None:
        assert self.g.nodes_areas(0) == 2.2
        assert self.g.nodes_areas(2) == 2.2
        npt.assert_equal(self.g.nodes_areas(), np.full(self.g.size, 2.2))

    def test_neighbors_count(self) -> None:
        assert self.g.neighbors_count(0) == 1
        assert self.g.neighbors_count(5) == 2

    def test_neighbors_indices(self) -> None:
        npt.assert_equal(self.g.neighbors_indices(0), np.array([1]))
        npt.assert_equal(self.g.neighbors_indices(5), np.array([4, 6]))

        with pytest.raises(IndexError, match="grid index out of range"):
            self.g.neighbors(11)

    def test_neighbors_distances(self) -> None:
        npt.assert_equal(self.g.neighbors_distances(0), np.array([2.2]))
        npt.assert_equal(self.g.neighbors_distances(5), np.array([2.2, 2.2]))

        with pytest.raises(IndexError, match="grid index out of range"):
            self.g.neighbors(11)

    def test_neighbors(self) -> None:
        assert self.g.neighbors(0) == [Neighbor(1, 2.2, NodeStatus.CORE)]
        assert self.g.neighbors(1) == [
            Neighbor(0, 2.2, NodeStatus.FIXED_VALUE),
            Neighbor(2, 2.2, NodeStatus.CORE),
        ]
        assert self.g.neighbors(6) == [
            Neighbor(5, 2.2, NodeStatus.FIXED_VALUE),
            Neighbor(7, 2.2, NodeStatus.CORE),
        ]

        with pytest.raises(IndexError, match="grid index out of range"):
            self.g.neighbors(11)
