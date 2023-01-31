import numpy as np
import pytest
from fastscapelib.algo import (
    compute_basins,
    compute_donors,
    compute_drainage_area_grid_d,
    compute_drainage_area_mesh_d,
    compute_receivers_d8_d,
    compute_stack,
    fill_sinks_flat,
    fill_sinks_sloped,
    find_pits,
)
from fastscapelib.grid import (
    NodeStatus,
    ProfileBoundaryStatus,
    ProfileGrid,
    RasterBoundaryStatus,
    RasterGrid,
)


@pytest.fixture(scope="session")
def braun_example():
    """Small drainage network example from Braun and Willet (2013) used
    for testing.

    """
    network = {}
    network["receivers"] = np.array([1, 4, 1, 6, 4, 4, 5, 4, 6, 7], dtype=np.intp)
    network["ndonors"] = np.array([0, 2, 0, 0, 3, 1, 2, 1, 0, 0], dtype=np.intp)
    network["donors"] = np.array(
        [
            [-1, -1, -1, -1, -1, -1, -1, -1],
            [0, 2, -1, -1, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1, -1, -1],
            [1, 5, 7, -1, -1, -1, -1, -1],
            [6, -1, -1, -1, -1, -1, -1, -1],
            [3, 8, -1, -1, -1, -1, -1, -1],
            [9, -1, -1, -1, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1, -1, -1],
        ],
        dtype=np.intp,
    )
    network["stack"] = np.array([4, 1, 0, 2, 5, 6, 3, 8, 7, 9])

    return network


def test_compute_receivers_d8():
    elevation = np.array(
        [
            [0.82, 0.16, 0.14, 0.20],
            [0.71, 0.97, 0.41, 0.09],
            [0.49, 0.01, 0.19, 0.38],
            [0.29, 0.82, 0.09, 0.88],
        ],
        dtype="d",
    )

    active_nodes = np.array(
        [
            [False, False, False, False],
            [False, True, True, False],
            [False, True, True, False],
            [False, False, False, False],
        ],
        dtype="bool",
    )

    receivers = np.ones((16), dtype=np.intp) * -1
    expected_receivers = np.array(
        [0, 1, 2, 3, 4, 9, 7, 7, 8, 9, 9, 11, 12, 13, 14, 15], dtype=np.intp
    )

    dist2receivers = np.ones((16), dtype="d") * -1
    expected_dist2receivers = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        dtype="d",
    )
    compute_receivers_d8_d(receivers, dist2receivers, elevation, active_nodes, 1.0, 1.0)

    np.testing.assert_array_equal(receivers, expected_receivers)
    np.testing.assert_allclose(dist2receivers, expected_dist2receivers)


def test_compute_donors(braun_example):
    ndonors = np.empty(10, dtype=np.intp)
    donors = np.ones((10, 8), dtype=np.intp) * -1

    compute_donors(ndonors, donors, braun_example["receivers"])

    np.testing.assert_array_equal(ndonors, braun_example["ndonors"])
    np.testing.assert_array_equal(donors, braun_example["donors"])


def test_compute_stack(braun_example):
    stack = np.empty(10, dtype=np.intp)
    compute_stack(
        stack,
        braun_example["ndonors"],
        braun_example["donors"],
        braun_example["receivers"],
    )

    np.testing.assert_array_equal(stack, braun_example["stack"])


def test_compute_basins(braun_example):
    basins = np.empty(10, dtype=np.intp)
    outlets_or_pits = np.ones(10, dtype=np.intp) * -1

    expected_basins = np.zeros_like(basins)
    expected_outlets_or_pits = np.array(
        [4, -1, -1, -1, -1, -1, -1, -1, -1, -1], dtype=np.intp
    )

    nbasins = compute_basins(
        basins, outlets_or_pits, braun_example["stack"], braun_example["receivers"]
    )

    assert nbasins == 1
    np.testing.assert_array_equal(basins, expected_basins)
    np.testing.assert_array_equal(outlets_or_pits, expected_outlets_or_pits)


def test_find_pits():
    # simple 4x4 test case with fixed boundaries and pit located at (2, 1)
    # 12 boundary nodes (i.e., 1-node open basin) + 1 pit = 13 basins
    nbasins = 13
    outlets_or_pits = np.array(
        [0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15, 9, -1, -1, -1], dtype=np.intp
    )
    active_nodes = np.array(
        [
            [False, False, False, False],
            [False, True, True, False],
            [False, True, True, False],
            [False, False, False, False],
        ],
        dtype="bool",
    )
    pits = np.ones(16, dtype=np.intp) * -1

    expected_pits = np.array(
        [9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1], dtype=np.intp
    )

    npits = find_pits(pits, outlets_or_pits, active_nodes, nbasins)

    assert npits == 1
    np.testing.assert_array_equal(pits, expected_pits)


def test_compute_drainage_area_mesh(braun_example):
    drainage_area = np.empty(10, dtype="d")
    cell_area = np.ones(10, dtype="d") * 2

    expected_drainage_area = np.array(
        [2.0, 6.0, 2.0, 2.0, 20.0, 8.0, 6.0, 4.0, 2.0, 2.0], dtype="d"
    )

    compute_drainage_area_mesh_d(
        drainage_area, cell_area, braun_example["stack"], braun_example["receivers"]
    )

    np.testing.assert_allclose(drainage_area, expected_drainage_area)


def test_compute_drainage_area_grid(braun_example):
    drainage_area = np.empty((2, 5), dtype="d")
    expected_drainage_area = np.array(
        [[2.0, 6.0, 2.0, 2.0, 20.0], [8.0, 6.0, 4.0, 2.0, 2.0]], dtype="d"
    )

    compute_drainage_area_grid_d(
        drainage_area, braun_example["stack"], braun_example["receivers"], 1.0, 2.0
    )

    np.testing.assert_allclose(drainage_area, expected_drainage_area)


@pytest.fixture(scope="session", params=["profile", "raster"])
def grid_type(request):
    return request.param


@pytest.fixture(scope="session")
def grid(grid_type):
    if grid_type == "profile":
        bs = ProfileBoundaryStatus(NodeStatus.FIXED_VALUE_BOUNDARY)
        return ProfileGrid(4, 1.0, bs, [])

    elif grid_type == "raster":
        bs = RasterBoundaryStatus(NodeStatus.FIXED_VALUE_BOUNDARY)
        return RasterGrid([3, 3], np.array([1.0, 1.0]), bs, [])


@pytest.fixture(scope="session")
def elevation(grid_type):
    if grid_type == "profile":
        return np.array([3.0, 0.1, 0.1, 2.0], dtype="d")

    elif grid_type == "raster":
        return np.array([[0.5, 0.4, 3.0], [3.0, 0.1, 3.0], [3.0, 3.0, 3.0]], dtype="d")


def test_fill_sinks_flat(grid_type, grid, elevation):
    res = elevation.copy()
    fill_sinks_flat(grid, res)

    if grid_type == "profile":
        expected = np.array([3.0, 2.0, 2.0, 2.0])
        np.testing.assert_array_equal(res, expected)

    elif grid_type == "raster":
        assert res[1, 1] == res[0, 1]
        assert res[0, 1] == elevation[0, 1]


def test_fill_sinks_sloped(grid_type, grid, elevation):
    res = elevation.copy()
    fill_sinks_sloped(grid, res)

    if grid_type == "profile":
        assert np.all(res[[0, 3]] == elevation[[0, 3]])
        assert res[2] > res[3]
        assert res[1] > res[2]

    elif grid_type == "raster":
        assert res[1, 1] >= res[0, 1]
        assert res[0, 1] == elevation[0, 1]
