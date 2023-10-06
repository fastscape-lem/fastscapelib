#include <cmath>
#include <limits>
#include <array>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "gtest/gtest.h"
#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xmath.hpp"

#include "fastscapelib/eroders/spl.hpp"
#include "fastscapelib/grid/raster_grid.hpp"
#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"
#include "fastscapelib/utils/utils.hpp"


using namespace xt::placeholders;
namespace fs = fastscapelib;


TEST(spl_eroder, ctor)
{
    using grid_type = fs::raster_grid;
    using size_type = grid_type::size_type;
    using shape_type = grid_type::shape_type;

    double spacing = 300;
    shape_type shape{ 2, 2 };

    auto grid = fs::raster_grid(shape, { spacing, spacing }, fs::node_status::fixed_value);

    auto flow_graph = fs::flow_graph<fs::raster_grid>(grid, { fs::single_flow_router() });

    double k_coef = 1e-3;
    double area_exp = 0.5;
    double slope_exp = 1.;
    xt::xtensor<double, 2> k_coef_arr = xt::ones<double>({ 2, 2 }) * k_coef;
    double tolerance = 1e-3;

    auto eroder = fs::make_spl_eroder(flow_graph, k_coef, area_exp, slope_exp, tolerance);

    EXPECT_TRUE(xt::all(xt::equal(eroder.k_coef(), xt::flatten(k_coef_arr))));
    EXPECT_EQ(eroder.area_exp(), area_exp);
    EXPECT_EQ(eroder.slope_exp(), slope_exp);
    EXPECT_EQ(eroder.tolerance(), tolerance);

    // setters
    eroder.set_k_coef(k_coef_arr);
    EXPECT_TRUE(xt::all(xt::equal(eroder.k_coef(), xt::flatten(k_coef_arr))));

    xt::xtensor<double, 2> k_coef_wrong_shape = xt::ones<double>({ 4, 5 });
    EXPECT_THROW(eroder.set_k_coef(k_coef_wrong_shape), std::runtime_error);

    eroder.set_area_exp(0.4);
    EXPECT_EQ(eroder.area_exp(), 0.4);

    eroder.set_slope_exp(1.2);
    EXPECT_EQ(eroder.slope_exp(), 1.2);

    {
        SCOPED_TRACE("error for n != 1 and multiple flow directions");

        auto flow_graph = fs::flow_graph<fs::raster_grid>(grid, { fs::multi_flow_router(1.0) });
        auto eroder = fs::make_spl_eroder(flow_graph, k_coef, area_exp, 1, tolerance);

        EXPECT_THROW(eroder.set_slope_exp(1.5), std::invalid_argument);
    }
}

/**
 * Setup for the evolution of a 1-d bedrock channel profile as a
 * result of uplift vs. erosion using the Stream-Power law.
 *
 * Drainage area along the profile is estimated using Hack's law.
 *
 */
class spl_eroder__profile_grid : public ::testing::Test
{
public:
    spl_eroder__profile_grid()
    {
        uplift(0) = 0.;
    }

protected:
    using flow_graph_type = fs::flow_graph<fs::profile_grid>;
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
    double area_exp = 0.5;
    double hack_coef = 6.69;
    double hack_exp = 1.67;

    using spl_eroder_type = fs::spl_eroder<flow_graph_type>;

    xt::xtensor<double, 1> drainage_area = hack_coef * xt::pow(x, hack_exp);

    // left base level
    fs::profile_boundary_status left_base_level{ fs::node_status::fixed_value,
                                                 fs::node_status::core };

    fs::profile_grid grid
        = fs::profile_grid(static_cast<size_type>(nnodes), spacing, left_base_level);

    template <class K>
    xt::xtensor<double, 1> get_steady_slope_numerical(flow_graph_type& flow_graph,
                                                      const K& k_coef,
                                                      double slope_exp,
                                                      double dt,
                                                      double tolerance);

    xt::xtensor<double, 1> get_steady_slope_analytical(double slope_exp);
};


/**
 * Get slope along the channel profile from the numerical solution at
 * steady state.
 */
template <class K>
xt::xtensor<double, 1>
spl_eroder__profile_grid::get_steady_slope_numerical(
    flow_graph_type& flow_graph, const K& k_coef, double slope_exp, double dt, double tolerance)
{
    auto eroder = spl_eroder_type(flow_graph, k_coef, area_exp, slope_exp, tolerance);

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
spl_eroder__profile_grid::get_steady_slope_analytical(double slope_exp)
{
    // exclude 1st node for proper comparison with the numerical solution
    auto drainage_area_ = xt::view(drainage_area, xt::range(1, _));

    xt::xtensor<double, 1> slope = (std::pow(u_rate / k_coef_scalar, 1. / slope_exp)
                                    * xt::pow(drainage_area_, -area_exp / slope_exp));

    return slope;
}


TEST_F(spl_eroder__profile_grid, update_routes)
{
    // not really testing SPL eroder but might useful to diagnose other tests
    auto flow_graph = flow_graph_type(grid, { fs::single_flow_router() });

    flow_graph.update_routes(elevation);

    xt::xtensor<index_t, 1> receivers = xt::arange<index_t>(-1, nnodes - 1);
    xt::xtensor<double, 1> dist2receivers = xt::ones<double>({ nnodes }) * spacing;
    xt::xtensor<index_t, 1> dfs_indices = xt::arange<index_t>(0, nnodes);
    receivers(0) = 0;
    dist2receivers(0) = 0.;

    EXPECT_TRUE(xt::all(xt::equal(receivers, xt::col(flow_graph.impl().receivers(), 0))));
    EXPECT_TRUE(
        xt::all(xt::equal(dist2receivers, xt::col(flow_graph.impl().receivers_distance(), 0))));
    EXPECT_TRUE(xt::all(xt::equal(dfs_indices, flow_graph.impl().dfs_indices())));
}


TEST_F(spl_eroder__profile_grid, erode__k_coef_scalar)
{
    auto flow_graph = flow_graph_type(grid, { fs::single_flow_router() });

    std::array<double, 3> slope_exp_vals{ 1., 2., 4. };

    for (const auto& slope_exp : slope_exp_vals)
    {
        SCOPED_TRACE("test against analytical solution of slope "
                     "along a 1-d profile at steady-state for slope_exp = "
                     + std::to_string(slope_exp));

        elevation = (length + x0 - x) * 1e-4;

        auto slope_n = get_steady_slope_numerical(flow_graph, k_coef_scalar, slope_exp, 1e4, 1e-3);
        auto slope_a = get_steady_slope_analytical(slope_exp);

        EXPECT_TRUE(xt::allclose(slope_n, slope_a, 1e-5, 1e-5));
        EXPECT_TRUE(n_corr == 0);
    }
}


TEST_F(spl_eroder__profile_grid, erode__k_coef_array)
{
    auto flow_graph = flow_graph_type(grid, { fs::single_flow_router() });

    std::array<double, 3> slope_exp_vals{ 1., 2., 4. };

    xt::xtensor<double, 1> k_coeff = xt::full_like(elevation, k_coef_scalar);

    for (const auto& slope_exp : slope_exp_vals)
    {
        SCOPED_TRACE("test against analytical solution of slope "
                     "using a 1-d array for k_coef and slope_exp = "
                     + std::to_string(slope_exp));

        elevation = (length + x0 - x) * 1e-4;

        auto slope_n = get_steady_slope_numerical(flow_graph, k_coeff, slope_exp, 1e4, 1e-3);
        auto slope_a = get_steady_slope_analytical(slope_exp);

        EXPECT_TRUE(xt::allclose(slope_n, slope_a, 1e-5, 1e-5));
        EXPECT_TRUE(n_corr == 0);
    }
}

TEST_F(spl_eroder__profile_grid, numerical_stability)
{
    auto flow_graph = flow_graph_type(grid, { fs::single_flow_router() });

    // crazy large values
    double k_coef = 1e8;
    double area_exp = 0.5;
    double slope_exp = 1;
    double dt = 1e8;

    auto eroder = spl_eroder_type(flow_graph, k_coef, area_exp, slope_exp);

    // crazy elevation with big hole
    elevation = (length + x0 - x) * 1e-4 + 1e3;
    elevation(0) = 0.0;
    elevation(10) = -1000;

    flow_graph.update_routes(elevation);

    {
        SCOPED_TRACE("linear case (n = 1)");

        eroder.erode(elevation, drainage_area, dt);
        EXPECT_GT(eroder.n_corr(), 0);
    }
    {
        SCOPED_TRACE("non-linear case");

        eroder.set_slope_exp(0.8);
        eroder.erode(elevation, drainage_area, dt);
        EXPECT_GT(eroder.n_corr(), 0);
    }
}


TEST(spl_eroder__raster_grid, tiny_grid)
{
    namespace fs = fastscapelib;

    using grid_type = fs::raster_grid;
    using flow_graph_type = fs::flow_graph<fs::raster_grid>;
    using size_type = grid_type::size_type;
    using shape_type = grid_type::shape_type;

    double spacing = 300;
    shape_type shape{ 2, 2 };

    // top border base-level
    fs::node_status fixed = fs::node_status::fixed_value;
    fs::node_status core = fs::node_status::core;
    fs::raster_boundary_status top_base_level{ { core, core, fixed, core } };

    auto grid = fs::raster_grid(shape, { spacing, spacing }, top_base_level);
    auto flow_graph = flow_graph_type(grid, { fs::single_flow_router() });

    double k_coef = 1e-3;
    xt::xtensor<double, 2> k_coef_arr = xt::ones<double>({ 2, 2 }) * k_coef;
    double area_exp = 0.5;
    double slope_exp = 1.;
    double tolerance = 1e-3;
    double dy = 300.;

    auto eroder = fs::make_spl_eroder(flow_graph, k_coef, area_exp, slope_exp, tolerance);

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
        xt::equal(xt::xtensor<size_type, 1>({ 0, 2, 1, 3 }), flow_graph.impl().dfs_indices())));

    EXPECT_TRUE(xt::all(xt::equal(xt::xtensor<double, 2>({ { 2 * a, 2 * a }, { a, a } }),
                                  flow_graph.accumulate(1.))));

    // spl erosion
    double dt = 1.;  // use small time step (compare with explicit scheme)

    auto erosion = eroder.erode(elevation, drainage_area, dt);

    double slope = h / dy;
    double err = dt * k_coef * std::pow(a, area_exp) * std::pow(slope, slope_exp);
    xt::xtensor<double, 2> expected_erosion = { { 0., 0. }, { err, err } };

    EXPECT_TRUE(xt::allclose(erosion, expected_erosion, 1e-5, 1e-5));
    EXPECT_TRUE(eroder.n_corr() == 0);
}


/*
 * High-level test: SPL is convervative and doesn't create new closed
 * depressions in the topography.
 */
TEST(spl_eroder, convervation_of_mass)
{
    using grid_type = fs::raster_grid;
    using flow_graph_type = fs::flow_graph<grid_type>;
    using size_type = typename grid_type::size_type;

    auto core = fs::node_status::core;
    auto fixed = fs::node_status::fixed_value;
    fs::raster_boundary_status bottom_fixed({ core, core, core, fixed });
    double spacing = 300.0;
    double cell_area = spacing * spacing;
    grid_type grid({ 50, 40 }, { spacing, spacing }, bottom_fixed);

    flow_graph_type graph_single(grid, { fs::pflood_sink_resolver(), fs::single_flow_router() });
    flow_graph_type graph_multi(grid, { fs::pflood_sink_resolver(), fs::multi_flow_router(0.0) });

    auto test_func = [&](flow_graph_type& graph, double slope_exp, double dt)
    {
        auto eroder = fs::make_spl_eroder(graph, 1e-3, 0.5, slope_exp, 1e-3);

        // planar surface tilted along the y-axis with random perturbations
        xt::xarray<double> elevation
            = (xt::random::rand<double>(grid.shape())
               * xt::view(xt::arange(grid.shape()[0]) * 2, xt::all(), xt::newaxis()));

        xt::xarray<double> uplift_rate = elevation / 1e3;

        xt::xarray<double> filled_elevation(elevation.shape());
        xt::xarray<double> drainage_area(elevation.shape());
        xt::xarray<double> erosion(elevation.shape());
        xt::xarray<double> sediment_acc(elevation.shape());

        size_t n_steps = 20;
        std::array<size_t, 1> shape{ n_steps };
        xt::xtensor<double, 1> actual(shape);
        xt::xtensor<double, 1> expected(shape);

        for (size_t i = 0; i < n_steps; i++)
        {
            elevation += uplift_rate * dt;
            filled_elevation = graph.update_routes(elevation);
            graph.accumulate(drainage_area, 1.0);
            erosion = eroder.erode(filled_elevation, drainage_area, dt);
            graph.accumulate(sediment_acc, erosion);

            // The sum of eroded material accumulated at base level
            // should be equal to the sum of total eroded material
            actual(i) = xt::sum(xt::row(sediment_acc, -1))();
            expected(i) = xt::sum(erosion * cell_area)();
        }

        EXPECT_TRUE(xt::allclose(actual, expected));
    };

    {
        SCOPED_TRACE("single flow router, n=1, dt=1e3");
        test_func(graph_single, 1.0, 1e3);
    }
    {
        SCOPED_TRACE("single flow router, n=2, dt=1e3");
        test_func(graph_single, 2.0, 1e3);
    }
    {
        SCOPED_TRACE("single flow router, n=1, dt=5e4");
        test_func(graph_single, 1.0, 5e4);
    }
    {
        SCOPED_TRACE("multi flow router, n=1, dt=1e3");
        test_func(graph_multi, 1.0, 1e3);
    }
    {
        SCOPED_TRACE("multi flow router, n=1, dt=5e4");
        test_func(graph_multi, 1.0, 5e4);
    }
}
