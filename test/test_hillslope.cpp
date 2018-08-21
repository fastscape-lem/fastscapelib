#include <array>

#include "gtest/gtest.h"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xtensor.hpp"

#include "fastscapelib/hillslope.hpp"


namespace fs = fastscapelib;


const double pi = 3.141592653589793238462643383279502884;


/* Get the fundamental solution of 2-d diffusion at time ``t``.
 *
 * The fundamental solution of diffusion has the form of a Gaussian
 * and is derived from assuming a dirac delta distribution as initial
 * conditions and an infinite domain.
 *
 * ``k_coef`` is assumed constant but might be an array.
 */
template<class X, class Y, class K>
auto solve_diffusion_analytical(X&& x, Y&& y, K&& k_coef, double t)
{
    auto fact = 4 * k_coef * t;
    xt::xtensor<double, 2> res = (xt::exp(-(xt::pow<2>(x) + xt::pow<2>(y)) / fact) /
                                  (fact * pi));

    return res;
}


template<class E1, class E2>
double compute_l2_norm(E1&& e1, E2&& e2)
{
    return 1. / static_cast<double>(e1.size()) * xt::sum(xt::pow<2>(e2 - e1))();
}


TEST(hillslope, erode_linear_diffusion)
{
    // test against fundamental solution of 2-d diffusion for multiple
    // values of dt (use l2-norm to compare results)
    auto grid = xt::meshgrid(xt::linspace<double>(-20, 20, 101),
                             xt::linspace<double>(-20, 20, 51));
    xt::xtensor<double, 2> y = std::get<0>(grid);
    xt::xtensor<double, 2> x = std::get<1>(grid);
    auto dy = 0.4;
    auto dx = 0.8;

    double k_coef = 1e-3;
    auto k_coef_arr = xt::full_like(x, k_coef);

    double t0 = 2e3;

    // initial conditions = fundamental solution of diffusion at time t0
    auto elevation_init = solve_diffusion_analytical(x, y, k_coef, t0);
    auto erosion = xt::empty_like(elevation_init);

    std::array<double, 4> dt_values {1e3, 1e4, 1e5, 1e6};

    // these are chosen arbitrarily after manual inspection (small enough values)
    std::array<double, 4> l2_norm_thresholds {1e-9, 1e-6, 1e-5, 1e-5};

    for (size_t k=0; k<4; ++k)
    {
        auto dt = dt_values[k];

        SCOPED_TRACE("test with scalar (double) k_coef and dt = " + std::to_string(dt));

        fs::erode_linear_diffusion(erosion, elevation_init, k_coef, dt, dx, dy);
        auto elevation_numerical = elevation_init - erosion;
        auto elevation_analytical = solve_diffusion_analytical(x, y, k_coef, t0 + dt);

        auto l2_norm = compute_l2_norm(elevation_analytical, elevation_numerical);

        EXPECT_TRUE(l2_norm < l2_norm_thresholds[k]);
    }

    for (size_t k=0; k<4; ++k)
    {
        auto dt = dt_values[k];

        SCOPED_TRACE("test with array k_coef and dt = " + std::to_string(dt));

        fs::erode_linear_diffusion(erosion, elevation_init, k_coef_arr, dt, dx, dy);
        auto elevation_numerical = elevation_init - erosion;
        auto elevation_analytical = solve_diffusion_analytical(x, y, k_coef, t0 + dt);

        auto l2_norm = compute_l2_norm(elevation_analytical, elevation_numerical);

        EXPECT_TRUE(l2_norm < l2_norm_thresholds[k]);
    }

}
