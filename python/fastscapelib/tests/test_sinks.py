import pytest
import numpy as np

from _fastscapelib_py.grid import (
    NodeStatus,
    RasterBoundaryStatus,
    RasterGrid,
    ProfileBoundaryStatus,
    ProfileGrid,
)
from _fastscapelib_py.sinks import fill_sinks_flat, fill_sinks_sloped


@pytest.fixture(scope='session', params=['profile', 'raster'])
def grid_type(request):
    return request.param


@pytest.fixture(scope='session')
def grid(grid_type):
    if grid_type == 'profile':
        bs = ProfileBoundaryStatus(NodeStatus.FIXED_VALUE_BOUNDARY)
        return ProfileGrid(4, 1.0, bs, [])

    elif grid_type == 'raster':
        bs = RasterBoundaryStatus(NodeStatus.FIXED_VALUE_BOUNDARY)
        return RasterGrid([3, 3], np.array([1.0, 1.0]), bs, [])


@pytest.fixture(scope='session')
def elevation(grid_type):
    if grid_type == 'profile':
        return np.array([3.0, 0.1, 0.1, 2.0], dtype='d')

    elif grid_type == 'raster':
        return np.array([[0.5, 0.4, 3.0],
                         [3.0, 0.1, 3.0],
                         [3.0, 3.0, 3.0]], dtype='d')


def test_fill_sinks_flat(grid_type, grid, elevation):
    res = elevation.copy()
    fill_sinks_flat(grid, res)

    if grid_type == 'profile':
        expected = np.array([3.0, 2.0, 2.0, 2.0])
        np.testing.assert_array_equal(res, expected)

    elif grid_type == 'raster':
        assert res[1, 1] == res[0, 1]
        assert res[0, 1] == elevation[0, 1]


def test_fill_sinks_sloped(grid_type, grid, elevation):
    res = elevation.copy()
    fill_sinks_sloped(grid, res)

    if grid_type == 'profile':
        assert np.all(res[[0, 3]] == elevation[[0, 3]])
        assert res[2] > res[3]
        assert res[1] > res[2]

    elif grid_type == 'raster':
        assert res[1, 1] >= res[0, 1]
        assert res[0, 1] == elevation[0, 1]
