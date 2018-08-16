import pytest
import numpy as np

import fastscapelib


def _solve_diffusion_analytical(x, y, k_coef, t):
    fact = 4 * k_coef * t
    return np.exp(-(x * x + y * y) / fact) / (fact * np.pi)


def _compute_l2_norm(a1, a2):
    return 1. / a1.size * np.sum(a2**2 - a1**2)


@pytest.mark.parametrize("k_coef_type", ["constant", "variable"])
def test_erode_linear_diffusion(k_coef_type):
    x, y = np.meshgrid(np.linspace(-20, 20, 51),
                       np.linspace(-20, 20, 101))
    dy = 0.4
    dx = 0.8

    k_coef = 1e-3
    dt = 1e3

    t0 = 2e3

    elevation_init = _solve_diffusion_analytical(x, y, k_coef, t0)
    erosion = np.empty_like(elevation_init)

    elevation_analytical = _solve_diffusion_analytical(x, y, k_coef, t0 + dt)

    if k_coef_type == 'constant':
        func = fastscapelib.erode_linear_diffusion_d
        k = k_coef
    elif k_coef_type == 'variable':
        func = fastscapelib.erode_linear_diffusion_var_d
        k = np.full_like(x, k_coef)

    func(erosion, elevation_init, k, dt, dx, dy)
    elevation_numerical = elevation_init - erosion

    l2_norm = _compute_l2_norm(elevation_analytical, elevation_numerical)

    assert l2_norm < 1e-9
