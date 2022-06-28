#include <cmath>
#include <limits>
#include <array>
#include <string>
#include <type_traits>

#include "gtest/gtest.h"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xio.hpp"

#include "fastscapelib/eroders/spl.hpp"
#include "fastscapelib/grid/raster_grid.hpp"
#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"
#include "fastscapelib/utils/utils.hpp"


using namespace xt::placeholders;
namespace fs = fastscapelib;


/**
 * Setup for the evolution of a 1-d bedrock channel profile as a
 * result of uplift vs. erosion using the Stream-Power law.
 *
 * Drainage area along the profile is estimated using Hack's law.
 *
 */
class erode_power_profile_grid : public ::testing::Test
{
public:
    erode_power_profile_grid()
    {
        uplift(0) = 0.;
    }

protected:
    using flow_graph_type
        = fs::flow_graph<fs::profile_grid, fs::single_flow_router, fs::no_sink_resolver>;
    using size_type = typename flow_graph_type::size_type;

    size_type n_corr = 0;
    int nnodes = 101;
    double spacing = 300.;
    double x0 = 300.;
    double length = (nnodes - 1) * spacing;

    xt::xtensor<double, 1> x
        = xt::linspace<double>(length + x0, x0, static_cast<std::size_t>(nnodes));
    xt::xtensor<double, 1> elevation = (length + x0 - x) * 1e-4;

    double u_rate = 1e-3;
    xt::xtensor<double, 1> uplift = xt::ones<double>({ nnodes }) * u_rate;

    double k_coef_scalar = 1e-3;
    double m_exp = 0.5;
    double hack_coef = 6.69;
    double hack_exp = 1.67;

    using spl_eroder_type = fs::spl_eroder<flow_graph_type>;

    xt::xtensor<double, 1> drainage_area = hack_coef * xt::pow(x, hack_exp);

    fs::profile_grid grid = fs::profile_grid(
        static_cast<size_type>(nnodes), spacing, fs::node_status::fixed_value_boundary);

    template <class K>
    xt::xtensor<double, 1> get_steady_slope_numerical(
        flow_graph_type& flow_graph, const K& k_coef, double n_exp, double dt, double tolerance);

    xt::xtensor<double, 1> get_steady_slope_analytical(double n_exp);
};


/**
 * Get slope along the channel profile from the numerical solution at
 * steady state.
 */
template <class K>
xt::xtensor<double, 1>
erode_power_profile_grid::get_steady_slope_numerical(
    flow_graph_type& flow_graph, const K& k_coef, double n_exp, double dt, double tolerance)
{
    auto eroder = spl_eroder_type(flow_graph, k_coef, m_exp, n_exp, tolerance);

    // run model until steady state is reached
    double elevation_diff = std::numeric_limits<double>::max();
    int iter = 0;

    while ((elevation_diff > 1e-3) && (iter < 1000))
    {
        xt::xtensor<double, 1> elevation_prev = elevation;
        elevation += uplift * dt;

        // routes should be invariant throughout the simulation
        // but better to update at each step in case there's a bug
        // that creates closed depressions.
        flow_graph.update_routes(elevation);

        auto erosion = eroder.erode(elevation, drainage_area, dt);

        n_corr += eroder.n_corr();

        elevation -= erosion;

        auto sum_squares = xt::sum(xt::pow(elevation - elevation_prev, 2))();
        elevation_diff = std::sqrt(sum_squares);
        ++iter;
    }

    // compute slope using backward finite difference scheme
    auto elevation_node = xt::view(elevation, xt::range(1, _));
    auto elevation_downstream = xt::view(elevation, xt::range(_, nnodes - 1));
    xt::xtensor<double, 1> slope = (elevation_node - elevation_downstream) / spacing;

    return slope;
}


/**
 * Get slope along the channel profile from the analytical solution at
 * steady state
 *
 * Note: this should be used to compare with the numerical solution only if
 * ``k_coef`` is invariant in space!
 */
xt::xtensor<double, 1>
erode_power_profile_grid::get_steady_slope_analytical(double n_exp)
{
    // exclude 1st node for proper comparison with the numerical solution
    auto drainage_area_ = xt::view(drainage_area, xt::range(1, _));

    xt::xtensor<double, 1> slope
        = (std::pow(u_rate / k_coef_scalar, 1. / n_exp) * xt::pow(drainage_area_, -m_exp / n_exp));

    return slope;
}


TEST_F(erode_power_profile_grid, flow_graph)
{
    // not really testing SPL eroder but might useful to diagnose other tests
    auto flow_graph = flow_graph_type(grid, fs::single_flow_router(), fs::no_sink_resolver());

    flow_graph.update_routes(elevation);

    xt::xtensor<index_t, 1> receivers = xt::arange<index_t>(-1, nnodes - 1);
    xt::xtensor<double, 1> dist2receivers = xt::ones<double>({ nnodes }) * spacing;
    xt::xtensor<index_t, 1> stack = xt::arange<index_t>(0, nnodes);
    receivers(0) = 0;
    dist2receivers(0) = 0.;

    EXPECT_TRUE(xt::all(xt::equal(receivers, xt::col(flow_graph.impl().receivers(), 0))));
    EXPECT_TRUE(
        xt::all(xt::equal(dist2receivers, xt::col(flow_graph.impl().receivers_distance(), 0))));
    EXPECT_TRUE(xt::all(xt::equal(stack, flow_graph.impl().dfs_stack())));
}


TEST_F(erode_power_profile_grid, scalar_k_coef)
{
    auto flow_graph = flow_graph_type(grid, fs::single_flow_router(), fs::no_sink_resolver());

    std::array<double, 3> n_exp_vals{ 1., 2., 4. };

    for (const auto& n_exp : n_exp_vals)
    {
        SCOPED_TRACE("test against analytical solution of slope "
                     "along a 1-d profile at steady-state for n_exp = "
                     + std::to_string(n_exp));

        elevation = (length + x0 - x) * 1e-4;

        auto slope_n = get_steady_slope_numerical(flow_graph, k_coef_scalar, n_exp, 1e4, 1e-3);
        auto slope_a = get_steady_slope_analytical(n_exp);

        EXPECT_TRUE(xt::allclose(slope_n, slope_a, 1e-5, 1e-5));
        EXPECT_TRUE(n_corr == 0);
    }
}


TEST_F(erode_power_profile_grid, array_k_coef)
{
    auto flow_graph = flow_graph_type(grid, fs::single_flow_router(), fs::no_sink_resolver());

    std::array<double, 3> n_exp_vals{ 1., 2., 4. };

    xt::xtensor<double, 1> k_coeff = xt::full_like(elevation, k_coef_scalar);

    for (const auto& n_exp : n_exp_vals)
    {
        SCOPED_TRACE("test against analytical solution of slope "
                     "using a 1-d array for k_coef and n_exp = "
                     + std::to_string(n_exp));

        elevation = (length + x0 - x) * 1e-4;

        auto slope_n = get_steady_slope_numerical(flow_graph, k_coeff, n_exp, 1e4, 1e-3);
        auto slope_a = get_steady_slope_analytical(n_exp);

        EXPECT_TRUE(xt::allclose(slope_n, slope_a, 1e-5, 1e-5));
        EXPECT_TRUE(n_corr == 0);
    }
}

TEST_F(erode_power_profile_grid, stability_threshold)
{
    auto flow_graph = flow_graph_type(grid, fs::single_flow_router(), fs::no_sink_resolver());

    {
        SCOPED_TRACE("test arbitrary limitation of erosion");

        elevation = (length + x0 - x) * 1e-4;

        // n_exp value super large -> massive (unstable) erosion
        auto slope_n = get_steady_slope_numerical(flow_graph, k_coef_scalar, 0.8, 1e4, 1e-3);
        EXPECT_TRUE(n_corr > 0);
    }
}


TEST(erode_stream_pow, raster_grid)
{
    namespace fs = fastscapelib;

    using grid_type = fs::raster_grid;
    using size_type = grid_type::size_type;
    using shape_type = grid_type::shape_type;

    {
        SCOPED_TRACE("test on a tiny (2x2) 2-d square grid "
                     "with a planar surface tilted in y (rows) "
                     "and set k_coef with a 2-d array");

        double spacing = 300;
        shape_type shape{ 2, 2 };

        auto grid
            = fs::raster_grid(shape, { spacing, spacing }, fs::node_status::fixed_value_boundary);

        auto flow_graph
            = fs::make_flow_graph(grid, fs::single_flow_router(), fs::no_sink_resolver());

        double k_coef = 1e-3;
        xt::xtensor<double, 2> k_coef_arr = xt::ones<double>({ 2, 2 }) * k_coef;
        double m_exp = 0.5;
        double n_exp = 1.;
        double tolerance = 1e-3;
        double dy = 300.;

        auto eroder = fs::make_spl_eroder(flow_graph, k_coef, m_exp, n_exp, tolerance);

        double h = 1.;
        double a = spacing * spacing;

        xt::xtensor<double, 2> elevation{ { 0., 0. }, { h, h } };

        // flow routing (tests useful to diagnose spl tests)
        flow_graph.update_routes(elevation);
        auto drainage_area = flow_graph.accumulate(1.);

        EXPECT_TRUE(xt::all(xt::equal(xt::xtensor<size_type, 1>({ 0, 1, 0, 1 }),
                                      xt::col(flow_graph.impl().receivers(), 0))));

        EXPECT_TRUE(xt::all(xt::equal(xt::xtensor<double, 1>({ 0., 0., spacing, spacing }),
                                      xt::col(flow_graph.impl().receivers_distance(), 0))));

        EXPECT_TRUE(xt::all(
            xt::equal(xt::xtensor<size_type, 1>({ 0, 2, 1, 3 }), flow_graph.impl().dfs_stack())));

        EXPECT_TRUE(xt::all(xt::equal(xt::xtensor<double, 2>({ { 2 * a, 2 * a }, { a, a } }),
                                      flow_graph.accumulate(xt::ones_like(elevation)))));

        // spl erosion
        double dt = 1.;  // use small time step (compare with explicit scheme)

        auto erosion = eroder.erode(elevation, drainage_area, dt);

        double slope = h / dy;
        double err = dt * k_coef * std::pow(a, m_exp) * std::pow(slope, n_exp);
        xt::xtensor<double, 2> expected_erosion = { { 0., 0. }, { err, err } };

        EXPECT_TRUE(xt::allclose(erosion, expected_erosion, 1e-5, 1e-5));
        EXPECT_TRUE(eroder.n_corr() == 0);
    }
}
