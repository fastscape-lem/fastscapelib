import numpy as np
import numpy.testing as npt
import pytest

from fastscapelib.grid import (
    Neighbor,
    NodeStatus,
    RasterBoundaryStatus,
    RasterGrid,
    RasterNeighbor,
)


class TestRasterBoundaryStatus:
    def setup_method(self) -> None:
        self.bs1 = RasterBoundaryStatus(
            [
                NodeStatus.CORE,
                NodeStatus.FIXED_VALUE,
                NodeStatus.CORE,
                NodeStatus.FIXED_VALUE,
            ]
        )
        self.bs2 = RasterBoundaryStatus(NodeStatus.CORE)
        self.bs3 = RasterBoundaryStatus(NodeStatus.LOOPED)

    def test__init__(self) -> None:
        # call the setup_method
        pass

    def test_left(self) -> None:
        assert self.bs1.left == NodeStatus.CORE
        assert self.bs2.left == NodeStatus.CORE
        assert self.bs3.left == NodeStatus.LOOPED

    def test_right(self) -> None:
        assert self.bs1.right == NodeStatus.FIXED_VALUE
        assert self.bs2.right == NodeStatus.CORE
        assert self.bs3.right == NodeStatus.LOOPED

    def test_bottom(self) -> None:
        assert self.bs1.bottom == NodeStatus.FIXED_VALUE
        assert self.bs2.bottom == NodeStatus.CORE
        assert self.bs3.bottom == NodeStatus.LOOPED

    def test_top(self) -> None:
        assert self.bs1.top == NodeStatus.CORE
        assert self.bs2.top == NodeStatus.CORE
        assert self.bs3.top == NodeStatus.LOOPED

    def test_is_horizontal_looped(self) -> None:
        assert not self.bs1.is_horizontal_looped
        assert not self.bs2.is_horizontal_looped
        assert self.bs3.is_horizontal_looped

    def test_is_vertical_looped(self) -> None:
        assert not self.bs1.is_vertical_looped
        assert not self.bs2.is_vertical_looped
        assert self.bs3.is_vertical_looped


class TestRasterNeighbor:
    def setup_method(self) -> None:
        self.n = RasterNeighbor(5, 0, 5, 1.35, NodeStatus.CORE)

    def test__init__(self) -> None:
        # call the setup_method
        assert self.n == RasterNeighbor(5, 0, 5, 1.35, NodeStatus.CORE)
        pass

    def test_flatten_idx(self) -> None:
        assert self.n.flatten_idx == 5

        self.n.flatten_idx = 8
        assert self.n.flatten_idx == 8

    def test_row(self) -> None:
        assert self.n.row == 0

        self.n.row = 3
        assert self.n.row == 3

    def test_col(self) -> None:
        assert self.n.col == 5

        self.n.col = 8
        assert self.n.col == 8

    def test_distance(self) -> None:
        assert self.n.distance == 1.35

        self.n.distance = 1.38
        assert self.n.distance == 1.38

    def test_status(self) -> None:
        assert self.n.status == NodeStatus.CORE

        self.n.status = NodeStatus.FIXED_VALUE
        assert self.n.status == NodeStatus.FIXED_VALUE


class TestRasterGrid:
    def setup_method(self) -> None:
        self.bs = bs = RasterBoundaryStatus(NodeStatus.FIXED_VALUE)
        self.g = RasterGrid([5, 10], [2.2, 2.4], bs, {(0, 5): NodeStatus.FIXED_VALUE})

    def test_static_props(self) -> None:
        assert RasterGrid.is_structured is True
        assert RasterGrid.is_uniform is True
        assert RasterGrid.n_neighbors_max == 8

    def test_constructor(self) -> None:
        g1 = RasterGrid([10, 10], [2.3, 2.1], self.bs, {(0, 5): NodeStatus.FIXED_VALUE})
        assert g1.size == 100
        npt.assert_almost_equal(g1.spacing, [2.3, 2.1])
        npt.assert_almost_equal(g1.length, [20.7, 18.9])

        with pytest.raises(IndexError):
            RasterGrid(
                [5, 10],
                [2.2, 2.4],
                self.bs,
                {(20, 255): NodeStatus.FIXED_VALUE},
            )

    @pytest.mark.parametrize(
        "bounds_status",
        [
            NodeStatus.FIXED_VALUE,
            [NodeStatus.FIXED_VALUE] * 4,
            RasterBoundaryStatus(NodeStatus.FIXED_VALUE),
        ],
    )
    def test_constructor_bounds_status(self, bounds_status) -> None:
        g = RasterGrid([2, 2], [2.0, 2.0], bounds_status)
        expected = np.ones((2, 2)) * NodeStatus.FIXED_VALUE.value
        npt.assert_array_equal(g.nodes_status(), expected)

    def test_from_length(self) -> None:
        g = RasterGrid.from_length(
            [11, 11], np.r_[23, 21], self.bs, {(0, 5): NodeStatus.FIXED_VALUE}
        )
        assert g.size == 121
        npt.assert_almost_equal(g.spacing, [2.3, 2.1])
        npt.assert_almost_equal(g.length, [23.0, 21.0])

    @pytest.mark.parametrize(
        "bounds_status",
        [
            NodeStatus.FIXED_VALUE,
            [NodeStatus.FIXED_VALUE] * 4,
            RasterBoundaryStatus(NodeStatus.FIXED_VALUE),
        ],
    )
    def test_from_length_bounds_status(self, bounds_status) -> None:
        g = RasterGrid.from_length([2, 2], [20.0, 20.0], bounds_status)
        expected = np.ones((2, 2)) * NodeStatus.FIXED_VALUE.value
        npt.assert_array_equal(g.nodes_status(), expected)

    def test_size(self) -> None:
        assert self.g.size == 50

    def test_shape(self) -> None:
        npt.assert_equal(self.g.shape, [5, 10])

    def test_spacing(self) -> None:
        npt.assert_equal(self.g.spacing, [2.2, 2.4])

    def test_nodes_indices(self) -> None:
        grid = RasterGrid([3, 3], [2.0, 2.0], self.bs)
        npt.assert_equal(grid.nodes_indices(), np.arange(grid.size))
        npt.assert_equal(
            grid.nodes_indices(NodeStatus.FIXED_VALUE), [0, 1, 2, 3, 5, 6, 7, 8]
        )
        npt.assert_equal(grid.nodes_indices(NodeStatus.CORE), [4])
        assert not len(grid.nodes_indices(NodeStatus.FIXED_GRADIENT))

    def test_nodes_status(self) -> None:
        grid = RasterGrid([3, 3], [2.0, 2.0], self.bs)

        npt.assert_equal(grid.nodes_status(), [[1, 1, 1], [1, 0, 1], [1, 1, 1]])
        # should return a copy
        assert not np.shares_memory(grid.nodes_status(), grid.nodes_status())

        assert grid.nodes_status(0) == NodeStatus.FIXED_VALUE
        assert grid.nodes_status(4) == NodeStatus.CORE

    def test_nodes_areas(self) -> None:
        area = 2.2 * 2.4
        assert self.g.nodes_areas(0) == area
        assert self.g.nodes_areas(10) == area
        npt.assert_equal(self.g.nodes_areas(), np.full(self.g.shape, area))

    def test_neighbors_count(self) -> None:
        assert self.g.neighbors_count(0) == 3
        assert self.g.neighbors_count(15) == 8

    def test_neighbors_indices(self) -> None:
        npt.assert_equal(self.g.neighbors_indices(0), np.array([1, 10, 11]))
        npt.assert_equal(
            self.g.neighbors_indices(15), np.array([4, 5, 6, 14, 16, 24, 25, 26])
        )

        with pytest.raises(IndexError, match="grid index out of range"):
            self.g.neighbors(51)

    def test_neighbors_raster_indices(self) -> None:
        assert self.g.neighbors_indices(0, 0) == [(0, 1), (1, 0), (1, 1)]
        assert self.g.neighbors_indices(1, 5) == [
            (0, 4),
            (0, 5),
            (0, 6),
            (1, 4),
            (1, 6),
            (2, 4),
            (2, 5),
            (2, 6),
        ]

        with pytest.raises(IndexError, match="grid index out of range"):
            self.g.neighbors(10, 10)

    def test_neighbors_distances(self) -> None:
        dist_diag = np.sqrt(2.2**2 + 2.4**2)
        npt.assert_equal(self.g.neighbors_distances(0), np.array([2.4, 2.2, dist_diag]))
        npt.assert_equal(
            self.g.neighbors_distances(15),
            np.array([dist_diag, 2.2, dist_diag, 2.4, 2.4, dist_diag, 2.2, dist_diag]),
        )

        with pytest.raises(IndexError, match="grid index out of range"):
            self.g.neighbors(51)

    def test_neighbors(self) -> None:
        dist_diag = np.sqrt(2.2**2 + 2.4**2)
        assert self.g.neighbors(0) == [
            Neighbor(1, 2.4, NodeStatus.FIXED_VALUE),
            Neighbor(10, 2.2, NodeStatus.FIXED_VALUE),
            Neighbor(11, dist_diag, NodeStatus.CORE),
        ]

        with pytest.raises(IndexError, match="grid index out of range"):
            self.g.neighbors(51)

    def test_raster_neighbors(self) -> None:
        dist_diag = np.sqrt(2.2**2 + 2.4**2)

        assert self.g.neighbors(0, 0) == [
            RasterNeighbor(1, 0, 1, 2.4, NodeStatus.FIXED_VALUE),
            RasterNeighbor(10, 1, 0, 2.2, NodeStatus.FIXED_VALUE),
            RasterNeighbor(11, 1, 1, dist_diag, NodeStatus.CORE),
        ]

        with pytest.raises(IndexError, match="grid index out of range"):
            self.g.neighbors(10, 10)
