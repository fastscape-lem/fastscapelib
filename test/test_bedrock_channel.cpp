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
 * Setup for the evolution of a 1-d bedrock channel profile as a
 * result of uplift vs. erosion using the Stream-Power law.
 *
 * Drainage area along the profile is estimated using Hack's law.
 */
struct ChannelProfile1D
{
    double k_coef = 1e-3;
    double m_exp = 0.5;
    double n_exp = 1;

    double u_rate = 1e-3;

    int nnodes = 101;
    double spacing = 300.;
    double x0 = 300.;

    double hack_coef = 6.69;
    double hack_exp = 1.67;

    index_t n_corr = 0;

    ChannelProfile1D()
    {
        // fixed boundary (outlet)
        receivers(0) = 0;
        dist2receivers(0) = 0.;
        uplift(0) = 0.;
    }

    xt::xtensor<double, 1> get_steady_slope_numerical(double dt, double tolerance);
    xt::xtensor<double, 1> get_steady_slope_analytical(void);

    private:
    double length = (nnodes - 1) * spacing;
    xt::xtensor<double, 1> x = xt::linspace<double>(length + x0, x0,
                                                    static_cast<size_t>(nnodes));

    xt::xtensor<double, 1> drainage_area = hack_coef * xt::pow(x, hack_exp);
    xt::xtensor<double, 1> elevation = (x - (length + x0)) * 0.001;

    xt::xtensor<index_t, 1> receivers = xt::arange<index_t>(-1, nnodes - 1);
    xt::xtensor<double, 1> dist2receivers = xt::ones<double>({nnodes}) * spacing;
    xt::xtensor<index_t, 1> stack = xt::arange<index_t>(0, nnodes);

    xt::xtensor<double, 1> uplift = xt::ones<double>({nnodes}) * u_rate;
    xt::xtensor<double, 1> erosion = xt::zeros<double>({nnodes});
};


/**
 * Get slope along the channel profile from the numerical solution at
 * steady state.
 */
xt::xtensor<double, 1> ChannelProfile1D::get_steady_slope_numerical(double dt,
                                                                    double tolerance)
{
    // run model until steady state is reached
    double elevation_diff = std::numeric_limits<double>::max();

    while (elevation_diff > 1e-3)
    {
        xt::xtensor<double, 1> elevation_prev = elevation;
        elevation += uplift * dt;

        n_corr += fs::erode_stream_power(erosion, elevation, stack,
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
 * Get slope along the channel profile from the analytical solution at
 * steady state.
 */
xt::xtensor<double, 1> ChannelProfile1D::get_steady_slope_analytical(void)
{
    // exclude 1st node for proper comparison with the numerical solution
    auto drainage_area_ = xt::view(drainage_area, xt::range(1, _));

    xt::xtensor<double, 1> slope = (std::pow(u_rate / k_coef, 1. / n_exp) *
                                    xt::pow(drainage_area_, -m_exp / n_exp));

    return slope;
}


TEST(bedrock_channel, erode_stream_power)
{
    std::array<double, 3> n_exp_vals {1., 2., 4.};

    for (const auto& n_exp : n_exp_vals)
    {
        SCOPED_TRACE("test against analytical solution of slope "
                     "along a 1-d profile at steady-state for n_exp = " +
                     std::to_string(n_exp));

        auto profile_1d = ChannelProfile1D();
        profile_1d.n_exp = n_exp;

        auto slope_n = profile_1d.get_steady_slope_numerical(1e4, 1e-3);
        auto slope_a = profile_1d.get_steady_slope_analytical();

        EXPECT_TRUE(xt::allclose(slope_n, slope_a, 1e-5, 1e-5));
        EXPECT_TRUE(profile_1d.n_corr == 0);
    }

    {
        SCOPED_TRACE("test arbitrary limitation of erosion");

        auto profile_1d = ChannelProfile1D();
        profile_1d.n_exp = 0.8;

        auto slope_n = profile_1d.get_steady_slope_numerical(1e4, 1e-3);

        EXPECT_TRUE(profile_1d.n_corr > 0);
    }

    {
        SCOPED_TRACE("test on a tiny (2x2) 2-d square grid "
                     "with a planar surface tilted in y (rows)");

        double k_coef = 1e-3;
        double m_exp = 0.5;
        double n_exp = 1.;

        double dy = 300.;

        xt::xtensor<index_t, 1> receivers {0, 1, 0, 1};
        xt::xtensor<double, 1> dist2receivers = {0., 0., dy, dy};
        xt::xtensor<index_t, 1> stack {0, 2, 1, 3};

        double a = dy * dy;
        xt::xtensor<double, 2> drainage_area {{2*a, 2*a}, {a, a}};

        double h = 1.;
        xt::xtensor<double, 2> elevation {{0., 0.}, {h, h}};

        xt::xtensor<double, 2> erosion = xt::empty<double>({2, 2});

        double dt = 1.;  // use small time step (compare with explicit scheme)
        double tolerance = 1e-3;

        index_t n_corr = fs::erode_stream_power(
            erosion, elevation, stack,
            receivers, dist2receivers, drainage_area,
            k_coef, m_exp, n_exp, dt, tolerance);

        double slope = h / dy;
        double err = dt * k_coef * std::pow(a, m_exp) * std::pow(slope, n_exp);
        xt::xtensor<double, 2> expected_erosion = {{0., 0.}, {err, err}};

        EXPECT_TRUE(xt::allclose(erosion, expected_erosion, 1e-5, 1e-5));
    }
}
