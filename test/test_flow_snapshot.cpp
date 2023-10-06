#include "xtensor/xtensor.hpp"
#include "xtensor/xstrided_view.hpp"

#include "gtest/gtest.h"

#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/flow_snapshot.hpp"
#include "fastscapelib/grid/raster_grid.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace testing
    {

        class flow_snapshot : public ::testing::Test
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

            // for the tests here we don't really care about the actual elevation values
            xt::xtensor<double, 2> elevation{ { 0.6, 0.6, 0.6, 0.6 },
                                              { 0.4, 0.4, 0.4, 0.4 },
                                              { 0.2, 0.2, 0.2, 0.2 },
                                              { 0.1, 0.0, 0.1, 0.1 } };
        };

        TEST_F(flow_snapshot, ctor)
        {
            auto snapshot = fs::flow_snapshot("test");
            EXPECT_EQ(snapshot.name(), "flow_snapshot");
            EXPECT_EQ(snapshot.snapshot_name(), "test");
            EXPECT_EQ(snapshot.save_graph(), true);
            EXPECT_EQ(snapshot.save_elevation(), false);
        }

        TEST_F(flow_snapshot, error)
        {
            // no flow routing operator defined before snapshot
            EXPECT_THROW(
                flow_graph_type(grid, { fs::flow_snapshot("a"), fs::single_flow_router() }),
                std::invalid_argument);
        }

        TEST_F(flow_snapshot, items)
        {
            auto graph = flow_graph_type(grid,
                                         {
                                             fs::pflood_sink_resolver(),
                                             fs::flow_snapshot("a", false, true),
                                             fs::single_flow_router(),
                                             fs::flow_snapshot("b", true, false),
                                         });

            ASSERT_EQ(graph.graph_snapshot_keys(), std::vector<std::string>{ "b" });
            ASSERT_EQ(graph.elevation_snapshot_keys(), std::vector<std::string>{ "a" });

            ASSERT_EQ(graph.graph_snapshot("b").size(), grid.size());
            ASSERT_EQ(graph.elevation_snapshot("a").size(), grid.size());
        }

        TEST_F(flow_snapshot, graph)
        {
            auto graph = flow_graph_type(grid,
                                         {
                                             fs::single_flow_router(),
                                             fs::flow_snapshot("a"),
                                             fs::mst_sink_resolver(),
                                             fs::flow_snapshot("b"),
                                         });

            graph.update_routes(elevation);
            auto& snapshot_a = graph.graph_snapshot("a");
            auto& snapshot_b = graph.graph_snapshot("b");

            // no operator in snapshot graphs
            ASSERT_EQ(snapshot_a.operators().size(), 0);
            ASSERT_EQ(snapshot_b.operators().size(), 0);

            ASSERT_NE(&graph, &snapshot_a);
            ASSERT_NE(&graph, &snapshot_b);

            // high-level interface test
            auto actual = graph.accumulate(1.0);
            auto expected = snapshot_b.accumulate(1.0);
            EXPECT_EQ(actual, expected);

            // test snapshot graphs are read-only
            ASSERT_THROW(snapshot_a.update_routes(elevation), std::runtime_error);
            ASSERT_THROW(snapshot_a.set_base_levels(std::vector<size_type>({ 0 })),
                         std::runtime_error);
            ASSERT_THROW(snapshot_a.set_mask(xt::zeros<bool>(grid.shape())), std::runtime_error);
        }

        TEST_F(flow_snapshot, elevation)
        {
            auto graph = flow_graph_type(grid,
                                         {
                                             fs::flow_snapshot("a", false, true),
                                             fs::pflood_sink_resolver(),
                                             fs::flow_snapshot("b", false, true),
                                             fs::single_flow_router(),
                                         });

            const auto& filled_elevation = graph.update_routes(elevation);
            const auto& snapshot_a = graph.elevation_snapshot("a");
            const auto& snapshot_b = graph.elevation_snapshot("b");

            EXPECT_EQ(elevation, snapshot_a);
            EXPECT_EQ(filled_elevation, snapshot_b);
        }
    }

}
