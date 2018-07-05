#include <cmath>
#include <limits>
#include <array>
#include <string>

#include "gtest/gtest.h"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xmath.hpp"

#include "fastscapelib/bedrock_channel.hpp"
#include "fastscapelib/utils.hpp"


using namespace xt::placeholders;
namespace fs = fastscapelib;


/**
 * Get slope of 1-d bedrock channel profile at steady state from
 * numerical simulation of uplift vs erosion using the stream-power
 * law.
 *
 * Drainage area is estimated using Hack's law.
 */
auto get_steady_slope_numerical(double k_coef,
                                double m_exp,
                                double n_exp,
                                double u_rate,
                                double hack_coef,
                                double hack_exp,
                                int nnodes,
                                double spacing,
                                double x0,
                                double dt,
                                double tolerance) -> xt::xtensor<double, 1>
{
    double length = (nnodes - 1) * spacing;
    auto x = xt::linspace<double>(length + x0, x0, static_cast<size_t>(nnodes));

    xt::xtensor<double, 1> drainage_area = hack_coef * xt::pow(x, hack_exp);
    xt::xtensor<double, 1> elevation = (x - (length + x0)) * 0.001;

    xt::xtensor<index_t, 1> receivers = xt::arange<index_t>(-1, nnodes - 1);
    xt::xtensor<double, 1> dist2receivers = xt::ones<double>({nnodes}) * spacing;
    xt::xtensor<index_t, 1> stack = xt::arange<index_t>(0, nnodes);

    xt::xtensor<double, 1> uplift = xt::ones<double>({nnodes}) * u_rate;
    xt::xtensor<double, 1> erosion = xt::zeros<double>({nnodes});

    // fixed boundary (outlet)
    receivers(0) = 0;
    dist2receivers(0) = 0.;
    uplift(0) = 0.;

    // run model until steady state is reached
    double elevation_diff = std::numeric_limits<double>::max();

    while (elevation_diff > 1e-3)
    {
        xt::xtensor<double, 1> elevation_prev = elevation;
        elevation += uplift * dt;

        index_t n_corr = fs::erode_stream_power(
            erosion, elevation, stack,
            receivers, dist2receivers, drainage_area,
            k_coef, m_exp, n_exp, dt, tolerance);

        elevation -= erosion;

        auto sum_squares = xt::sum(xt::pow(elevation - elevation_prev, 2))();
        elevation_diff = std::sqrt(sum_squares);
    }

    // compute slope using backward finite difference scheme
    auto elevation_node = xt::view(elevation, xt::range(1, _));
    auto elevation_downstream = xt::view(elevation, xt::range(_, nnodes - 1));
    xt::xtensor<double, 1> slope = (elevation_node - elevation_downstream) / spacing;

    return slope;
}


/**
 * Get slope of 1-d bedrock channel profile at steady state from
 * analytical solution of uplift vs. erosion using the stream-power
 * law.
 *
 * Drainage area is estimated using Hack's law.
 */
auto get_steady_slope_analytical(double k_coef,
                                 double m_exp,
                                 double n_exp,
                                 double u_rate,
                                 double hack_coef,
                                 double hack_exp,
                                 int nnodes,
                                 double spacing,
                                 double x0) -> xt::xtensor<double, 1>
{
    double length = (nnodes - 1) * spacing;
    auto x = xt::linspace<double>(length + x0, x0, static_cast<size_t>(nnodes));

    auto drainage_area = hack_coef * xt::pow(x, hack_exp);

    // exclude 1st node for comparison with numerical solution
    auto drainage_area_ = xt::view(drainage_area, xt::range(1, _));

    xt::xtensor<double, 1> slope = (std::pow(u_rate / k_coef, 1. / n_exp) *
                                    xt::pow(drainage_area_, -m_exp / n_exp));

    return slope;
}


TEST(bedrock_channel, erode_stream_power)
{
    // Test numerical vs. analytical solution of the slope of a 1-d
    // river profile at steady state.
    // Test for multiple n_exp values
    // TODO: add test for n_exp < 1 (not working currently)
    double k_coef = 1e-3;
    double m_exp = 0.5;
    std::array<double, 3> n_exp_vals {1., 2., 4.};

    double u_rate = 0.001;

    int nnodes = 101;
    double spacing = 300.;
    double x0 = 300.;

    double hack_coef = 6.69;
    double hack_exp = 1.67;

    double dt = 1e4;
    double tolerance = 1e-3;

    for (const auto& n_exp : n_exp_vals)
    {
        SCOPED_TRACE("steady-state analytical 1-d with n=" + std::to_string(n_exp));

        auto slope_n = get_steady_slope_numerical(k_coef, m_exp, n_exp,
                                                  u_rate, hack_coef, hack_exp,
                                                  nnodes, spacing, x0,
                                                  dt, tolerance);

        auto slope_a = get_steady_slope_analytical(k_coef, m_exp, n_exp,
                                                   u_rate, hack_coef, hack_exp,
                                                   nnodes, spacing, x0);

        EXPECT_TRUE(xt::allclose(slope_n, slope_a, 1e-5, 1e-5));
    }

    // TODO: Test arbitrary limitation of erosion

    // Test on a tiny (2x2) 2-d square grid with a planar surface
    // tilted in y (rows) and with all outlets on the 1st row.
    {
        SCOPED_TRACE("simple 2-d test");

        xt::xtensor<index_t, 1> receivers {0, 1, 0, 1};
        xt::xtensor<double, 1> dist2receivers = {0., 0., spacing, spacing};
        xt::xtensor<index_t, 1> stack {0, 2, 1, 3};

        double a = spacing * spacing;
        xt::xtensor<double, 2> drainage_area {{2*a, 2*a}, {a, a}};

        double h = 1.;
        xt::xtensor<double, 2> elevation {{0., 0.}, {h, h}};

        xt::xtensor<double, 2> erosion = xt::empty<double>({2, 2});

        double n_exp = 1.;
        double dt = 1.;  // use small time step (compare with explicit scheme)

        index_t n_corr = fs::erode_stream_power(
            erosion, elevation, stack,
            receivers, dist2receivers, drainage_area,
            k_coef, m_exp, n_exp, dt, tolerance);

        double slope = h / spacing;
        double err = dt * k_coef * std::pow(a, m_exp) * std::pow(slope, n_exp);
        xt::xtensor<double, 2> expected_erosion = {{0., 0.}, {err, err}};

        EXPECT_TRUE(xt::allclose(erosion, expected_erosion, 1e-5, 1e-5));
    }
}
