#include <cmath>
#include <limits>
#include <unordered_set>
#include <vector>

#include "gtest/gtest.h"
#include "xtensor/xtensor.hpp"

#include "fastscapelib/flow/basin_graph.hpp"
#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/grid/raster_grid.hpp"

#include "gtest/gtest.h"


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace testing
    {

        class basin_graph : public ::testing::TestWithParam<fs::mst_method>
        {
        protected:
            using grid_type = fs::raster_grid;
            using flow_graph_type = fs::flow_graph<grid_type>;
            using flow_graph_impl_type = flow_graph_type::impl_type;
            using basin_graph_type = fs::basin_graph<flow_graph_impl_type>;

            using size_type = typename basin_graph_type::size_type;
            using data_type = typename basin_graph_type::data_type;

            using edge_type = typename basin_graph_type::edge;
            using edges_type = std::vector<edge_type>;

            // bottom border base-level
            fs::node_status fixed = fs::node_status::fixed_value;
            fs::node_status core = fs::node_status::core;
            fs::raster_boundary_status bottom_base_level{ { core, core, core, fixed } };

            grid_type grid = grid_type({ 5, 5 }, { 1.0, 1.0 }, bottom_base_level);

            // planar surface tilted along the y-axis
            // + one closed depression with pit at row/col index (2, 2)
            //   and pass between (2, 2) and (3, 1)
            xt::xtensor<double, 2> elevation{ { 1.00, 1.00, 1.00, 1.00, 1.00 },
                                              { 0.70, 0.70, 0.70, 0.20, 0.70 },
                                              { 0.50, 0.50, 0.10, 0.50, 0.50 },
                                              { 0.20, 0.19, 0.20, 0.20, 0.20 },
                                              { 0.00, 0.00, 0.00, 0.00, 0.00 } };

            flow_graph_type fgraph = flow_graph_type(grid, { fs::single_flow_router() });

            void setup()
            {
                fgraph.update_routes(elevation);
                // triggers computing the basins
                fgraph.basins();
            }

            void expect_edges_equal(edges_type& left, edges_type& right)
            {
                EXPECT_EQ(left.size(), right.size());

                for (size_t i = 0; i < left.size(); ++i)
                {
                    EXPECT_TRUE(left[i] == right[i]) << "edges are not equal at index " << i;
                }
            }
        };


        TEST_P(basin_graph, ctor)
        {
            setup();

            basin_graph_type bgraph = basin_graph_type(fgraph.impl(), GetParam());

            EXPECT_EQ(bgraph.perf_boruvka(), 0);

            {
                SCOPED_TRACE("test outlets");

                EXPECT_EQ(bgraph.basins_count(), fgraph.impl().outlets().size());

                auto actual = bgraph.outlets();
                auto expected = fgraph.impl().outlets();
                EXPECT_TRUE(
                    std::equal(actual.begin(), actual.end(), expected.begin(), expected.end()));
            }
        }


        TEST_F(basin_graph, compute_basins)
        {
            // not really part of basin graph and already tested elsewhere but
            // useful for readability of the tests here and as quick diagnostic
            // in case any test below fails

            fgraph.update_routes(elevation);

            auto actual = fgraph.basins();
            xt::xtensor<size_type, 2> expected{ { 1, 0, 0, 0, 0 },
                                                { 1, 0, 0, 0, 0 },
                                                { 1, 0, 0, 0, 5 },
                                                { 1, 2, 3, 4, 5 },
                                                { 1, 2, 3, 4, 5 } };

            EXPECT_TRUE(xt::all(xt::equal(actual, expected)));
        }


        TEST_F(basin_graph, connect_basins)
        {
            setup();

            basin_graph_type bgraph = basin_graph_type(fgraph.impl(), fs::mst_method::kruskal);

            bgraph.update_routes(elevation);

            // raster grid diagonal spacing
            data_type diag = std::sqrt(2);

            // for edges with no assigned pass
            size_type init_idx = static_cast<size_type>(-1);
            data_type init_elev = std::numeric_limits<data_type>::lowest();

            auto actual = bgraph.edges();
            edges_type expected = { // all edges from inner basin 0 to outer basins
                                    { { 2, 0 }, { 16, 12 }, 0.19, diag },
                                    { { 0, 3 }, { 12, 17 }, 0.20, 1.0 },
                                    { { 0, 4 }, { 12, 18 }, 0.20, diag },
                                    { { 0, 1 }, { 11, 10 }, 0.50, 1.0 },
                                    { { 0, 5 }, { 8, 14 }, 0.50, diag },
                                    // all outer basins connected to the virtual
                                    // root basin (assigned to basin 1 as it is
                                    // the 1st parsed outer basin)
                                    { { 1, 2 }, { init_idx, init_idx }, init_elev, 0.0 },
                                    { { 1, 3 }, { init_idx, init_idx }, init_elev, 0.0 },
                                    { { 1, 4 }, { init_idx, init_idx }, init_elev, 0.0 },
                                    { { 1, 5 }, { init_idx, init_idx }, init_elev, 0.0 }
            };

            expect_edges_equal(actual, expected);
        }

        TEST_F(basin_graph, connect_basins_mask)
        {
            // masking all nodes of basin 5
            xt::xtensor<bool, 2> mask{ { false, false, false, false, false },
                                       { false, false, false, false, false },
                                       { false, false, false, false, true },
                                       { false, false, false, false, true },
                                       { false, false, false, false, true } };

            setup();
            fgraph.set_mask(mask);
            fgraph.update_routes(elevation);

            basin_graph_type bgraph = basin_graph_type(fgraph.impl(), fs::mst_method::kruskal);
            bgraph.update_routes(elevation);

            // copied from connect_basins test, removed edges with basin 5
            data_type diag = std::sqrt(2);
            size_type init_idx = static_cast<size_type>(-1);
            data_type init_elev = std::numeric_limits<data_type>::lowest();
            auto actual = bgraph.edges();
            edges_type expected = { // all edges from inner basin 0 to outer basins
                                    { { 2, 0 }, { 16, 12 }, 0.19, diag },
                                    { { 0, 3 }, { 12, 17 }, 0.20, 1.0 },
                                    { { 0, 4 }, { 12, 18 }, 0.20, diag },
                                    { { 0, 1 }, { 11, 10 }, 0.50, 1.0 },
                                    // all outer basins connected to the virtual
                                    // root basin (assigned to basin 1 as it is
                                    // the 1st parsed outer basin)
                                    { { 1, 2 }, { init_idx, init_idx }, init_elev, 0.0 },
                                    { { 1, 3 }, { init_idx, init_idx }, init_elev, 0.0 },
                                    { { 1, 4 }, { init_idx, init_idx }, init_elev, 0.0 }
            };

            expect_edges_equal(actual, expected);
        }

        TEST_P(basin_graph, compute_mst)
        {
            setup();

            basin_graph_type bgraph = basin_graph_type(fgraph.impl(), GetParam());

            bgraph.update_routes(elevation);

            // only the 1st edge shown in the previous test, i.e., with link (2,
            // 0), is kept in the tree since it has the lowest pass for inner
            // basin 0.
            //
            // all outer basins <-> virtual root basin edges (i.e., the 4 last
            // edges shown in the previous test) are also kept in the tree.
            //
            // the MST method (kruskal vs. boruvka) may not yield edges positions
            // in the same order -> compare unordered collections
            //
            using set_type = std::unordered_set<size_type>;

            const auto tree = bgraph.tree();
            set_type actual(tree.begin(), tree.end());
            set_type expected{ 0, 5, 6, 7, 8 };

            EXPECT_EQ(actual, expected);
        }


        TEST_F(basin_graph, orient_edges)
        {
            setup();

            basin_graph_type bgraph = basin_graph_type(fgraph.impl(), fs::mst_method::kruskal);

            bgraph.update_routes(elevation);

            // change the orientation of the (2, 0) edge (kept in the MST tree)
            // and re-run orient_edges
            edge_type e = bgraph.m_edges[0];

            std::swap(e.link[0], e.link[1]);
            std::swap(e.pass[0], e.pass[1]);

            bgraph.orient_edges();

            // check that the edge is correctly re-oriented
            //
            // (2, 0) is the counter flow direction, i.e., from (outer) basin 2
            // to (inner) basin 0
            edge_type actual = bgraph.edges()[0];
            edge_type expected{ { 2, 0 }, { 16, 12 }, 0.19, std::sqrt(2) };
            EXPECT_TRUE(actual == expected);
        }


        INSTANTIATE_TEST_SUITE_P(mst_methods,
                                 basin_graph,
                                 ::testing::Values(fs::mst_method::kruskal,
                                                   fs::mst_method::boruvka));

    }
}
