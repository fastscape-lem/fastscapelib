#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"
#include "fastscapelib/grid/raster_grid.hpp"

#include "xtensor/xtensor.hpp"
#include "xtensor/xio.hpp"

#include "gtest/gtest.h"


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace testing
    {

        class flow_graph : public ::testing::Test
        {
        protected:
            using flow_graph_type
                = fs::flow_graph<fs::raster_grid, fs::single_flow_router, fs::no_sink_resolver>;
            using grid_type = fs::raster_grid;
            using size_type = typename grid_type::size_type;

            // bottom border base-level
            fs::node_status fixed = fs::node_status::fixed_value_boundary;
            fs::node_status core = fs::node_status::core;
            fs::raster_boundary_status bottom_base_level{ { core, core, core, fixed } };

            grid_type grid = grid_type({ 4, 4 }, { 1.0, 1.0 }, bottom_base_level);

            // planar surface tilted along the y-axis + small carved channel
            xt::xtensor<double, 2> elevation{ { 0.6, 0.6, 0.6, 0.6 },
                                              { 0.4, 0.4, 0.4, 0.4 },
                                              { 0.2, 0.2, 0.2, 0.2 },
                                              { 0.1, 0.0, 0.1, 0.1 } };
        };

        TEST_F(flow_graph, ctor)
        {
            auto graph
                = fs::make_flow_graph(grid, fs::single_flow_router(), fs::no_sink_resolver());

            EXPECT_EQ(graph.size(), 16u);
        }

        TEST_F(flow_graph, update_routes)
        {
            auto graph
                = fs::make_flow_graph(grid, fs::single_flow_router(), fs::no_sink_resolver());
            graph.update_routes(elevation);
            graph.update_routes(elevation);  // check there is not memory effect

            xt::xtensor<size_type, 1> expected_fixed_receivers{ 4,  5,  6,  7,  8,  9,  10, 11,
                                                                13, 13, 13, 15, 12, 13, 14, 15 };

            EXPECT_TRUE(
                xt::all(xt::equal(xt::col(graph.impl().receivers(), 0), expected_fixed_receivers)));
        }

        TEST_F(flow_graph, accumulate)
        {
            auto graph
                = fs::make_flow_graph(grid, fs::single_flow_router(), fs::no_sink_resolver());

            graph.update_routes(elevation);

            xt::xarray<double> acc = xt::empty_like(elevation);
            xt::xtensor<double, 2> src = xt::ones_like(elevation);
            xt::xtensor<double, 2> expected{
                { 1., 1., 1., 1. }, { 2., 2., 2., 2. }, { 3., 3., 3., 3. }, { 1., 10., 1., 4. }
            };

            graph.accumulate(acc, 1.);
            EXPECT_TRUE(xt::allclose(expected, acc));

            graph.accumulate(acc, src);
            EXPECT_TRUE(xt::allclose(expected, acc));

            EXPECT_TRUE(xt::allclose(expected, graph.accumulate(src)));
            EXPECT_TRUE(xt::allclose(expected, graph.accumulate(1.)));
        }
    }
}
