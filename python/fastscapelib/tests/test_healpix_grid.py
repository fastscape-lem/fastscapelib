import numpy as np
import pytest

from fastscapelib.grid import HealpixGrid, NodeStatus

pytestmark = pytest.mark.skipif(
    HealpixGrid is None, reason="Fastscapelib not compiled with Healpix support"
)


@pytest.fixture
def nside() -> int:
    return 32


@pytest.fixture
def size() -> int:
    return 12288


@pytest.fixture
def nodes_status(size) -> np.ndarray:
    ghost = NodeStatus.GHOST
    return np.zeros(size, dtype=np.uint8)


class TestHealpixGrid:
    def test_static_properties(self) -> None:
        assert HealpixGrid.is_structured is False
        assert HealpixGrid.is_uniform is False
        assert HealpixGrid.n_neighbors_max == 8

    def test_constructor(self, nside, nodes_status, size) -> None:
        grid = HealpixGrid(nside, nodes_status)
        assert grid.nside == nside
        assert grid.size == size
        assert grid.radius == 6.371e6

    def test_nodes_status(self, nside, nodes_status) -> None:
        grid = HealpixGrid(nside, nodes_status)

        np.testing.assert_array_equal(grid.nodes_status(), nodes_status)

        nodes_status2 = nodes_status.copy()
        nodes_status2[0] = NodeStatus.LOOPED
        with pytest.raises(ValueError, match="node_status::looped is not allowed"):
            HealpixGrid(nside, nodes_status2)

        with pytest.raises(ValueError, match="invalid shape for nodes_status array"):
            HealpixGrid(nside, nodes_status[:10])

    def test_nodes_areas(self, nside, nodes_status) -> None:
        grid = HealpixGrid(nside, nodes_status)

        expected = 41509152987.45021
        assert abs(grid.nodes_areas(0) - expected) < 1e-5

        expected_arr = np.ones(grid.size) * expected
        np.testing.assert_allclose(grid.nodes_areas(), expected_arr)

    def test_nodes_lonlat(self, nside, nodes_status) -> None:
        grid = HealpixGrid(nside, nodes_status)

        expected = (0.7853981633974483, 1.5452801164374776)
        assert grid.nodes_lonlat(0) == expected
        lon, lat = grid.nodes_lonlat()
        assert (lon[0], lat[0]) == expected

    def test_nodes_xyz(self, nside, nodes_status) -> None:
        grid = HealpixGrid(nside, nodes_status)

        expected = (114937.4753788243, 114937.4753788243, 6368926.106770833)
        assert grid.nodes_xyz(0) == expected
        x, y, z = grid.nodes_xyz()
        assert (x[0], y[0], z[0]) == expected

    def test_neighbors_count(self, nside, nodes_status) -> None:
        grid = HealpixGrid(nside, nodes_status)

        assert grid.neighbors_count(0) == 7
        assert grid.neighbors_count(1984) == 6

    def test_neighbors_indices(self, nside, nodes_status) -> None:
        grid = HealpixGrid(nside, nodes_status)

        np.testing.assert_array_equal(
            grid.neighbors_indices(0), [11, 3, 2, 1, 6, 5, 13]
        )

    def test_neighbors_distances(self, nside, nodes_status) -> None:
        grid = HealpixGrid(nside, nodes_status)

        expected = np.array(
            [
                0.04752047,
                0.03608146,
                0.05102688,
                0.03608146,
                0.04752047,
                0.02914453,
                0.0510435,
            ],
        )
        np.testing.assert_allclose(grid.neighbors_distances(0), expected, rtol=1e-5)
