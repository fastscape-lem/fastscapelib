import numpy as np

import fastscapelib


def test_erode_stream_power():
    # Test on a tiny (2x2) 2-d square grid with a planar surface
    # tilted in y (rows) and with all outlets on the 1st row.
    spacing = 300.

    receivers = np.array([0, 1, 0, 1], dtype=np.intp)
    dist2receivers = np.array([0., 0., spacing, spacing], dtype='d')
    stack = np.array([0, 2, 1, 3], dtype=np.intp)

    a = spacing**2
    drainage_area = np.array([[2 * a, 2 * a], [a, a]], dtype='d')

    h = 1.
    elevation = np.array([[0., 0.], [h, h]], dtype='d')

    erosion = np.empty_like(elevation)

    k_coef = 1e-3
    m_exp = 0.5
    n_exp = 1.

    dt = 1   # use small time step (compare with explicit scheme)
    tolerance = 1e-3

    n_corr = fastscapelib.erode_stream_power_d(
        erosion, elevation, stack,
        receivers, dist2receivers, drainage_area,
        k_coef, m_exp, n_exp, dt, tolerance)

    slope = h / spacing
    err = dt * k_coef * a**m_exp * slope**n_exp
    expected_erosion = np.array([[0., 0.], [err, err]], dtype='d')

    np.testing.assert_allclose(erosion, expected_erosion, atol=1e-5)
    assert n_corr == 0
