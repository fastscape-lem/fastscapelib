#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/grid/raster_grid.hpp"

#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

#include "gtest/gtest.h"


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace testing
    {

        class multi_flow_router : public ::testing::Test
        {
        protected:
            using size_type = typename fs::raster_grid::size_type;
            using grid_type = fs::raster_grid;
            using flow_graph_type = fs::flow_graph<grid_type>;

            // test on a 3x3 tiny grid with base levels at all border nodes
            grid_type grid = grid_type({ 3, 3 }, { 1.0, 1.0 }, fs::node_status::fixed_value);

            double sqrt2 = std::sqrt(2.0);

            xt::xtensor<double, 2> elevation{ { 0.0, sqrt2, 0.0 },
                                              { sqrt2, sqrt2 * 2, sqrt2 },
                                              { 0.0, sqrt2, 0.0 } };
        };

        TEST_F(multi_flow_router, ctor)
        {
            auto router = fs::multi_flow_router(1.0);
            EXPECT_EQ(router.name(), "multi_flow_router");
            EXPECT_EQ(router.m_slope_exp, 1.0);

            router.m_slope_exp = 0.0;
            EXPECT_EQ(router.m_slope_exp, 0.0);
        }

        TEST_F(multi_flow_router, graph_topology)
        {
            flow_graph_type graph(grid, { fs::multi_flow_router(0.0) });
            graph.update_routes(elevation);

            {
                SCOPED_TRACE("receivers count");
                auto actual = graph.impl().receivers_count();
                xt::xtensor<size_type, 1> expected{ 1, 1, 1, 1, 8, 1, 1, 1, 1 };
                EXPECT_EQ(actual, expected);
            }
            {
                SCOPED_TRACE("receivers");
                auto actual = xt::row(graph.impl().receivers(), 4);
                xt::xtensor<size_type, 1> expected{ 0, 1, 2, 3, 5, 6, 7, 8 };
                EXPECT_EQ(actual, expected);
            }
            {
                SCOPED_TRACE("receivers distance");
                auto actual = xt::row(graph.impl().receivers_distance(), 4);
                xt::xtensor<double, 1> expected{ sqrt2, 1.0, sqrt2, 1.0, 1.0, sqrt2, 1.0, sqrt2 };
                EXPECT_EQ(actual, expected);
            }
            {
                SCOPED_TRACE("donors count");
                auto actual = graph.impl().donors_count();
                xt::xtensor<size_type, 1> expected{ 1, 1, 1, 1, 0, 1, 1, 1, 1 };
                EXPECT_EQ(actual, expected);
            }
            {
                SCOPED_TRACE("donors");
                auto actual = xt::view(graph.impl().donors(), xt::keep(0, 1, 2, 3, 5, 6, 7, 8), 0);
                xt::xtensor<size_type, 1> expected{ 4, 4, 4, 4, 4, 4, 4, 4 };
                EXPECT_EQ(xt::flatten(actual), expected);
            }
        }

        TEST_F(multi_flow_router, mask)
        {
            flow_graph_type graph(grid, { fs::multi_flow_router(0.0) });

            {
                SCOPED_TRACE("masked node");

                // center node is masked, edge nodes are base levels
                xt::xtensor<bool, 2> mask{ { false, false, false },
                                           { false, true, false },
                                           { false, false, false } };
                graph.set_mask(mask);
                graph.update_routes(elevation);

                // no graph edge (each node has itself as receiver)
                auto actual = graph.impl().receivers_count();
                xt::xtensor<size_type, 1> expected{ 1, 1, 1, 1, 1, 1, 1, 1, 1 };
                EXPECT_EQ(actual, expected);
            }
            {
                SCOPED_TRACE("masked neighbor");

                // one masked edge node
                xt::xtensor<bool, 2> mask{ { false, false, false },
                                           { true, false, false },
                                           { false, false, false } };
                graph.set_mask(mask);
                graph.update_routes(elevation);

                // masked node not present in center node receivers
                auto actual = xt::row(graph.impl().receivers(), 4);
                xt::xtensor<size_type, 1> expected{ 0, 1, 2, 5, 6, 7, 8 };
                EXPECT_EQ(xt::view(actual, xt::range(0, 7)), expected);
            }
        }

        TEST_F(multi_flow_router, flow_partition)
        {
            {
                SCOPED_TRACE("slope_exp = 0");

                flow_graph_type graph(grid, { fs::multi_flow_router(0.0) });
                graph.update_routes(elevation);

                auto actual = xt::row(graph.impl().receivers_weight(), 4);
                xt::xtensor<double, 1> expected{ 1. / 8, 1. / 8, 1. / 8, 1. / 8,
                                                 1. / 8, 1. / 8, 1. / 8, 1. / 8 };
                EXPECT_TRUE(xt::allclose(actual, expected));
            }
            {
                SCOPED_TRACE("slope_exp = 2");

                flow_graph_type graph(grid, { fs::multi_flow_router(2.0) });
                graph.update_routes(elevation);

                auto actual = xt::row(graph.impl().receivers_weight(), 4);
                xt::xtensor<double, 1> expected{ 2. / 12, 1. / 12, 2. / 12, 1. / 12,
                                                 1. / 12, 2. / 12, 1. / 12, 2. / 12 };
                EXPECT_TRUE(xt::allclose(actual, expected));
            }
        }
    }
}
