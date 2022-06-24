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

            fs::node_status fb = fs::node_status::fixed_value_boundary;
            fs::raster_boundary_status fixed_value_status{ fb };

            grid_type grid = grid_type({ 4, 4 }, { 1.1, 1.2 }, fixed_value_status);

            xt::xtensor<double, 2> elevation{ { 0.82, 0.16, 0.14, 0.20 },
                                              { 0.71, 0.97, 0.41, 0.09 },
                                              { 0.49, 0.01, 0.19, 0.38 },
                                              { 0.29, 0.82, 0.09, 0.88 } };
        };

        TEST_F(flow_graph, ctor)
        {
            auto graph
                = fs::make_flow_graph(grid, fs::single_flow_router(), fs::no_sink_resolver());

            EXPECT_EQ(graph.grid().size(), 16u);  // dummy test
        }

        TEST_F(flow_graph, update_routes)
        {
            auto graph
                = fs::make_flow_graph(grid, fs::single_flow_router(), fs::no_sink_resolver());
            graph.update_routes(elevation);
            graph.update_routes(elevation);  // check there is not memory effect

            xt::xtensor<size_type, 1> expected_fixed_receivers{ 1, 2, 7, 7, 9, 9, 7, 7,
                                                                9, 9, 9, 7, 9, 9, 9, 14 };

            EXPECT_TRUE(
                xt::all(xt::equal(xt::col(graph.receivers(), 0), expected_fixed_receivers)));
        }

        TEST_F(flow_graph, accumulate)
        {
            auto graph
                = fs::make_flow_graph(grid, fs::single_flow_router(), fs::no_sink_resolver());
            graph.update_routes(elevation);

            xt::xtensor<double, 2> data1 = xt::ones_like(elevation);
            xt::xtensor<double, 2> data2{ { 1.1, 1.0, 1.1, 1.0 },
                                          { 1.1, 1.0, 1.1, 1.0 },
                                          { 1.1, 1.0, 1.1, 1.0 },
                                          { 1.1, 1.0, 1.1, 1.0 } };

            xt::xtensor<double, 2> expected1{ { 1.32, 2.64, 3.96, 1.32 },
                                              { 1.32, 1.32, 1.32, 9.24 },
                                              { 1.32, 11.88, 1.32, 1.32 },
                                              { 1.32, 1.32, 2.64, 1.32 } };

            xt::xtensor<double, 2> expected2{ { 1.452, 2.772, 4.224, 1.32 },
                                              { 1.452, 1.32, 1.452, 9.636 },
                                              { 1.452, 12.54, 1.452, 1.32 },
                                              { 1.452, 1.32, 2.772, 1.32 } };

            EXPECT_TRUE(xt::allclose(expected1, graph.accumulate(data1)));
            EXPECT_TRUE(xt::allclose(expected1, graph.accumulate(1.)));
            EXPECT_TRUE(xt::allclose(expected1 * 5., graph.accumulate(5.)));
            EXPECT_TRUE(xt::allclose(expected2, graph.accumulate(data2)));
        }
    }
}
