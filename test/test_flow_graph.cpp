#include <algorithm>

#include "xtensor/xtensor.hpp"
#include "xtensor/xstrided_view.hpp"

#include "gtest/gtest.h"

#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/flow_snapshot.hpp"
#include "fastscapelib/grid/raster_grid.hpp"


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace testing
    {

        class flow_graph : public ::testing::Test
        {
        protected:
            using flow_graph_type = fs::flow_graph<fs::raster_grid>;
            using grid_type = fs::raster_grid;
            using size_type = typename grid_type::size_type;

            // bottom border base-level
            fs::node_status fixed = fs::node_status::fixed_value;
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
            auto graph = flow_graph_type(grid, { fs::single_flow_router() });

            EXPECT_EQ(graph.size(), 16u);
            EXPECT_TRUE(xt::same_shape(graph.grid_shape(), grid.shape()));

            // no operator updating the graph (flow router)
            EXPECT_THROW(flow_graph_type(grid, { fs::flow_snapshot("s") }), std::invalid_argument);
        }

        TEST_F(flow_graph, base_levels)
        {
            auto graph = flow_graph_type(grid, { fs::single_flow_router() });

            {
                SCOPED_TRACE("default base levels");

                auto actual = graph.base_levels();
                // base_levels are not ordered (implementation detail)
                std::sort(actual.begin(), actual.end());
                EXPECT_EQ(actual, std::vector<size_type>({ 12, 13, 14, 15 }));
            }
            {
                SCOPED_TRACE("set base levels");

                graph.set_base_levels(std::vector<size_type>({ 5 }));
                EXPECT_EQ(graph.base_levels(), std::vector<size_type>({ 5 }));
            }
        }

        TEST_F(flow_graph, mask)
        {
            auto graph = flow_graph_type(grid, { fs::single_flow_router() });

            {
                SCOPED_TRACE("uninitialized mask (default)");

                EXPECT_EQ(graph.mask(), xt::xarray<bool>());
            }
            {
                SCOPED_TRACE("set mask");

                graph.set_mask(xt::ones<bool>(grid.shape()));
                EXPECT_EQ(graph.mask(), xt::ones<bool>(grid.shape()));
            }
            {
                SCOPED_TRACE("set mask error (shape mismatch)");

                EXPECT_THROW(graph.set_mask(xt::ones<bool>({ 10. })), std::runtime_error);
            }
        }

        TEST_F(flow_graph, operators)
        {
            auto graph
                = flow_graph_type(grid, { fs::single_flow_router(), fs::flow_snapshot("s") });

            std::vector<std::string> expected{ "single_flow_router", "flow_snapshot" };
            std::vector<std::string> actual;

            for (const auto* op : graph.operators())
            {
                actual.push_back(op->name());
            }

            EXPECT_EQ(actual, expected);
        }

        TEST_F(flow_graph, single_flow)
        {
            {
                SCOPED_TRACE("single flow is true");

                auto graph = flow_graph_type(grid, { fs::single_flow_router() });
                EXPECT_EQ(graph.single_flow(), true);
            }
            {
                SCOPED_TRACE("single flow is false");

                auto graph = flow_graph_type(grid, { fs::multi_flow_router(1.0) });
                EXPECT_EQ(graph.single_flow(), false);
            }
        }

        TEST_F(flow_graph, update_routes)
        {
            auto graph = flow_graph_type(grid, { fs::single_flow_router() });

            auto new_elevation = graph.update_routes(elevation);

            // re-run to check there is no memory effect
            graph.update_routes(elevation);

            EXPECT_TRUE(xt::all(xt::equal(new_elevation, elevation)));

            // more in-depth testing is done in flow router (methods) tests
            auto actual = xt::col(graph.impl().receivers(), 0);
            xt::xtensor<size_type, 1> expected{ 4,  5,  6,  7,  8,  9,  10, 11,
                                                13, 13, 13, 15, 12, 13, 14, 15 };

            EXPECT_EQ(actual, expected);
        }

        TEST_F(flow_graph, accumulate)
        {
            auto graph = flow_graph_type(grid, { fs::single_flow_router() });

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

        TEST_F(flow_graph, basins)
        {
            auto graph = flow_graph_type(grid, { fs::single_flow_router() });

            graph.update_routes(elevation);

            auto actual = graph.basins();

            xt::xtensor<size_type, 2> expected{
                { 1, 1, 1, 3 }, { 1, 1, 1, 3 }, { 1, 1, 1, 3 }, { 0, 1, 2, 3 }
            };

            EXPECT_EQ(actual, expected);

            EXPECT_TRUE(xt::all(xt::equal(xt::flatten(actual), graph.impl().basins())));

            {
                SCOPED_TRACE("with mask");

                xt::xtensor<bool, 2> mask{ { false, false, false, true },
                                           { false, false, false, true },
                                           { false, false, false, true },
                                           { false, false, false, true } };

                graph.set_mask(mask);
                graph.update_routes(elevation);

                auto actual = graph.basins();

                size_type no_basin = std::numeric_limits<size_type>::max();

                xt::xtensor<size_type, 2> expected{ { 1, 1, 1, no_basin },
                                                    { 1, 1, 1, no_basin },
                                                    { 1, 1, 1, no_basin },
                                                    { 0, 1, 2, no_basin } };

                EXPECT_EQ(actual, expected);
            }
        }
    }
}
