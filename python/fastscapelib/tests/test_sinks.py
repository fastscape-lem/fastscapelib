import numpy as np

import fastscapelib


def test_fill_sinks_flat():
    arr = np.array([[1.0, 2.0, 3.0],
                    [2.0, 0.1, 7.0],
                    [2.0, 5.0, 7.0]], dtype='d')

    fastscapelib.fill_sinks_flat_d(arr)

    assert arr[1, 1] == arr[0, 0]


def test_fill_sinks_sloped():
    arr = np.array([[1.0, 2.0, 3.0],
                    [2.0, 0.1, 7.0],
                    [2.0, 5.0, 7.0]], dtype='d')

    fastscapelib.fill_sinks_sloped_d(arr)

    assert arr[1, 1] > arr[0, 0]
