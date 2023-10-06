#include "xtensor/xbuilder.hpp"
#include "xtensor/xrandom.hpp"
#include "xtensor/xset_operation.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

#include "gtest/gtest.h"

#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"
#include "fastscapelib/grid/raster_grid.hpp"

namespace fs = fastscapelib;
using namespace xt::placeholders;


namespace fastscapelib
{

    // trivial router subclasses for parametrized tests
    struct multi_flow_router_0 : public fs::multi_flow_router
    {
        multi_flow_router_0()
            : fs::multi_flow_router(0.0)
        {
        }
    };

    struct multi_flow_router_2 : public fs::multi_flow_router
    {
        multi_flow_router_2()
            : fs::multi_flow_router(2.0)
        {
        }
    };

    namespace detail
    {
        template <class FG>
        class flow_operator_impl<FG, multi_flow_router_0, flow_graph_fixed_array_tag>
            : public flow_operator_impl<FG, multi_flow_router, flow_graph_fixed_array_tag>
        {
        public:
            using base_type = flow_operator_impl<FG, multi_flow_router, flow_graph_fixed_array_tag>;

            flow_operator_impl(std::shared_ptr<multi_flow_router_0> ptr)
                : base_type(std::move(ptr)){};
        };

        template <class FG>
        class flow_operator_impl<FG, multi_flow_router_2, flow_graph_fixed_array_tag>
            : public flow_operator_impl<FG, multi_flow_router, flow_graph_fixed_array_tag>
        {
        public:
            using base_type = flow_operator_impl<FG, multi_flow_router, flow_graph_fixed_array_tag>;

            flow_operator_impl(std::shared_ptr<multi_flow_router_2> ptr)
                : base_type(std::move(ptr)){};
        };
    }

    namespace testing
    {

        class flow_router_base : public ::testing::Test
        {
        protected:
            using grid_type = fs::raster_grid;
            using flow_graph_type = fs::flow_graph<grid_type>;
            using size_type = typename grid_type::size_type;

            // bottom border base-level
            fs::node_status fixed = fs::node_status::fixed_value;
            fs::node_status core = fs::node_status::core;
            fs::raster_boundary_status bottom_base_level{ { core, core, core, fixed } };

            fs::node_status fb = fs::node_status::fixed_value;
            fs::raster_boundary_status fixed_value_status{ fb };

            size_type nrows = 10;
            size_type ncols = 8;

            grid_type grid = grid_type({ nrows, ncols }, { 1.0, 1.0 }, bottom_base_level);

            // planar surface tilted along the y-axis with random perturbations
            xt::xarray<double> elevation
                = (xt::random::rand<double>(grid.shape())
                   * xt::view(xt::arange(nrows) * 2, xt::all(), xt::newaxis()));
        };

        template <class OP>
        class flow_router : public flow_router_base
        {
        protected:
            using router_type = OP;
        };

        using router_types
            = ::testing::Types<fs::single_flow_router, multi_flow_router_0, multi_flow_router_2>;
        TYPED_TEST_SUITE(flow_router, router_types);

        /*
         * High-level test: convervative flow routing
         */
        TYPED_TEST(flow_router, convervation_of_area)
        {
            using router_type = typename TestFixture::router_type;
            using flow_graph_type = typename TestFixture::flow_graph_type;

            // we want all the flow be routed to the last row (y-axis) -> resolve sinks
            auto graph = flow_graph_type(this->grid, { fs::pflood_sink_resolver(), router_type() });
            graph.update_routes(this->elevation);
            auto drainage_area = graph.accumulate(1.0);

            double actual = xt::sum(xt::row(drainage_area, -1))();

            // assumes grid cell (uniform) area is 1
            double expected = double(this->grid.size());

            EXPECT_NEAR(actual, expected, 1e-5);
        }

        /*
         * High-level test: DFS traversal and monotonic elevation
         */
        TYPED_TEST(flow_router, monotonic_dfs)
        {
            using router_type = typename TestFixture::router_type;
            using flow_graph_type = typename TestFixture::flow_graph_type;
            using size_type = typename TestFixture::size_type;

            // this test requires no closed depressions
            auto graph = flow_graph_type(this->grid, { fs::pflood_sink_resolver(), router_type() });
            const auto& filled_elevation = graph.update_routes(this->elevation);

            const auto& receivers_count = graph.impl().receivers_count();
            const auto& receivers = graph.impl().receivers();
            const auto& dfs_indices = graph.impl().dfs_indices();

            {
                SCOPED_TRACE("test forward iterator");

                // traverse the graph, check that elevation increases if
                // previous node is a direct receiver of current visited node
                // (dfs_indices always in bottom->up direction)
                size_type i = 0;
                for (const auto& inode : graph.impl().nodes_indices_bottomup())
                {
                    auto irec = receivers(inode, 0);
                    if (irec == inode)
                    {
                        i++;
                        continue;
                    }

                    auto iprev = dfs_indices(i - 1);
                    auto r_count = receivers_count(inode);

                    if (xt::any(xt::isin(iprev, xt::view(receivers, xt::range(_, r_count)))))
                    {
                        EXPECT_GE(filled_elevation.flat(inode), filled_elevation.flat(iprev));
                    }

                    i++;
                }
            }
            {
                SCOPED_TRACE("test reverse iterator");

                // traverse the graph in the top->down direction, check that
                // elevation decrease for the next node to visit unless it is an
                // outlet
                size_type i = this->grid.size();
                auto indices = graph.impl().nodes_indices_bottomup();

                for (auto inode = indices.rbegin(); inode != indices.rend(); ++inode)
                {
                    auto irec = receivers(*inode, 0);
                    if (irec == *inode)
                    {
                        i--;
                        continue;
                    }

                    auto inext = dfs_indices(i - 1);
                    auto irec_next = receivers(inext, 0);
                    if (inext != irec_next)
                    {
                        EXPECT_GE(filled_elevation.flat(*inode), filled_elevation.flat(inext));
                    }

                    i--;
                }
            }
        }
    }
}
