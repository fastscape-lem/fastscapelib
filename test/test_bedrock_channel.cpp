#include <limits>

#include "gtest/gtest.h"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xmath.hpp"

#include "fastscapelib/bedrock_channel.hpp"


using namespace xt::placeholders;
namespace fs = fastscapelib;


/**
 * Get slope of 1-d bedrock channel profile at steady state from
 * numerical simulation.
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

    // iterate until steady state is reached
    double elevation_diff = std::numeric_limits<double>::max();

    while (elevation_diff > 1e-3)
    {
        xt::xtensor<double, 1> elevation_prev = elevation;
        elevation += uplift * dt;
        fs::erode_stream_power(erosion, elevation, stack,
                               receivers, dist2receivers, drainage_area,
                               k_coef, m_exp, n_exp, dt, tolerance);
        elevation -= erosion;

        auto sum_squares = xt::sum(xt::pow(elevation - elevation_prev, 2))();
        elevation_diff = std::sqrt(sum_squares);
    }

    // compute slope
    auto elevation_right = xt::view(elevation, xt::range(1, _));
    auto elevation_left = xt::view(elevation, xt::range(_, nnodes - 1));
    xt::xtensor<double, 1> slope = (elevation_right - elevation_left) / spacing;

    return slope;
}


/**
 * Get slope of 1-d bedrock channel profile at steady state from
 * analytical solution.
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

    auto _slope = (std::pow(u_rate / k_coef, 1. / n_exp) *
                   xt::pow(drainage_area, -m_exp / n_exp));

    xt::xtensor<double, 1> slope = xt::view(_slope, xt::range(1, _));

    return slope;
}


TEST(bedrock_channel, erode_stream_power)
{
    double k_coef = 1e-3;
    double m_exp = 0.5;

    double u_rate = 0.001;

    int nnodes = 101;
    double spacing = 300.;
    double x0 = 300.;

    double hack_coef = 6.69;
    double hack_exp = 1.67;

    double dt = 1e4;
    double tolerance = 1e-3;

    {
        // Test against analytical solution of slope at steady state for
        // 1d river profile (n == 1)
        SCOPED_TRACE("steady-state analytical 1d n=1");

        double n_exp = 1.;

        auto slope_n = get_steady_slope_numerical(k_coef, m_exp, n_exp,
                                                  u_rate, hack_coef, hack_exp,
                                                  nnodes, spacing, x0,
                                                  dt, tolerance);

        auto slope_a = get_steady_slope_analytical(k_coef, m_exp, n_exp,
                                                   u_rate, hack_coef, hack_exp,
                                                   nnodes, spacing, x0);

        EXPECT_TRUE(xt::allclose(slope_n, slope_a, 1e-5, 1e-6));
    }

    // {
    //     // TODO: (regression test: fix bug)
    //     // Test against analytical solution of slope at steady state for
    //     // 1d river profile ( n != 1, check for Newton iterations)
    //     SCOPED_TRACE("steady-state analytical 1d n!=1");

    //     double n_exp = 1.5;

    //     auto slope_n = get_steady_slope_numerical(k_coef, m_exp, n_exp,
    //                                               u_rate, hack_coef, hack_exp,
    //                                               nnodes, spacing, x0,
    //                                               dt, tolerance);

    //     auto slope_a = get_steady_slope_analytical(k_coef, m_exp, n_exp,
    //                                                u_rate, hack_coef, hack_exp,
    //                                                nnodes, spacing, x0);

    //     EXPECT_TRUE(xt::allclose(slope_n, slope_a, 1e-5, 1e-6));
    // }

    // TODO: test for 2D (maybe just simple dumb test to ensure no
    //       error at compile/runtime is raised when using 2-d arrays)

    // Example in Braun and Willet, 2013 as a test case
    // - fixed drainage area = 2 for all cells
    // - fixed dist2receivers = 1 for all nodes (except outlet)
    /*xt::xtensor<index_t, 1> receivers {1, 4, 1, 6, 4, 4, 5, 4, 6, 7};
    xt::xtensor<double, 1> dist2receivers = xt::ones<double>({10}) * 2.;
    xt::xtensor<index_t, 1> stack {4, 1, 0, 2, 5, 6, 3, 8, 7, 9};
    xt::xtensor<double, 1> drainage_area = xt::ones<double>({10}) * 2.;
    xt::xtensor<double, 1> elevation = xt::zeros<double>({10});
    xt::xtensor<double, 1> erosion = xt::zeros<double>({10});

    fs::erode_stream_power(erosion, elevation, stack,
                           receivers, dist2receivers, drainage_area,
                           1e-7, 0.4, 1., 1000., 1e-3);
    */
    //EXPECT_EQ(arr(1, 1), arr(0, 0));
}
