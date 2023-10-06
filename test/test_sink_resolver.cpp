#include <cmath>
#include <limits>
#include <memory>

#include "xtensor/xarray.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xoperation.hpp"
#include "xtensor/xrandom.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

#include "gtest/gtest.h"

#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"
#include "fastscapelib/grid/raster_grid.hpp"


namespace fs = fastscapelib;

namespace fastscapelib
{

    namespace testing
    {

        class sink_resolver : public ::testing::Test
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

            grid_type grid = grid_type({ 5, 5 }, { 1.0, 1.0 }, bottom_base_level);

            // planar surface tilted along the y-axis
            // + one closed depression with pit at row/col index (1, 3)
            //   and pass between (2, 2) and (3, 1)
            xt::xarray<double> elevation{ { 1.00, 1.00, 1.00, 1.00, 1.00 },
                                          { 0.70, 0.70, 0.70, 0.10, 0.70 },
                                          { 0.50, 0.50, 0.15, 0.50, 0.50 },
                                          { 0.20, 0.19, 0.20, 0.20, 0.20 },
                                          { 0.00, 0.00, 0.00, 0.00, 0.00 } };
        };


        TEST_F(sink_resolver, no_sink_resolver)
        {
            auto graph = flow_graph_type(grid, { fs::single_flow_router() });
            const auto& new_elevation = graph.update_routes(elevation);

            {
                SCOPED_TRACE("test unchanged elevation");
                EXPECT_TRUE(xt::all(xt::equal(new_elevation, elevation)));
            }

            {
                SCOPED_TRACE("test flow still trapped in pit");
                EXPECT_EQ(graph.impl().receivers()(8, 0), 8);
            }
        }


        TEST_F(sink_resolver, pflood_sink_resolver)
        {
            EXPECT_EQ(fs::pflood_sink_resolver().name(), "pflood_sink_resolver");

            auto graph
                = flow_graph_type(grid, { fs::pflood_sink_resolver(), fs::single_flow_router() });
            const auto& new_elevation = graph.update_routes(elevation);

            {
                SCOPED_TRACE("test filled elevation (tiny slope)");

                xt::xtensor<bool, 2> mask_not_filled = xt::ones<bool>(elevation.shape());
                mask_not_filled(1, 3) = false;
                mask_not_filled(2, 2) = false;

                auto new_unchanged = xt::where(mask_not_filled, new_elevation, 0.0);
                auto unchanged = xt::where(mask_not_filled, elevation, 0.0);
                EXPECT_TRUE(xt::all(xt::equal(new_unchanged, unchanged)));

                EXPECT_TRUE(new_elevation(2, 2) > new_elevation(3, 1));
                EXPECT_TRUE(new_elevation(1, 3) > new_elevation(2, 2));
            }

            {
                SCOPED_TRACE("test resolved flow");

                EXPECT_EQ(graph.impl().receivers()(8, 0), 12);
                EXPECT_EQ(graph.impl().receivers_distance()(8, 0), std::sqrt(2.0));
                EXPECT_EQ(graph.impl().receivers()(12, 0), 16);
                EXPECT_EQ(graph.impl().receivers_distance()(12, 0), std::sqrt(2.0));
            }
        }


        class mst_sink_resolver
            : public sink_resolver
            , public ::testing::WithParamInterface<fs::mst_method>
        {
        };


        TEST_F(mst_sink_resolver, ctor)
        {
            {
                SCOPED_TRACE("test default constructor");

                auto resolver = fs::mst_sink_resolver();
                EXPECT_EQ(resolver.name(), "mst_sink_resolver");
                EXPECT_EQ(resolver.m_basin_method, fs::mst_method::kruskal);
                EXPECT_EQ(resolver.m_route_method, fs::mst_route_method::carve);
            }

            {
                SCOPED_TRACE("test parameters constructor");

                auto resolver
                    = fs::mst_sink_resolver(fs::mst_method::boruvka, fs::mst_route_method::carve);
                EXPECT_EQ(resolver.m_basin_method, fs::mst_method::boruvka);
                EXPECT_EQ(resolver.m_route_method, fs::mst_route_method::carve);
            }

            {
                SCOPED_TRACE("test initializer list");

                auto resolver
                    = fs::mst_sink_resolver{ fs::mst_method::boruvka, fs::mst_route_method::basic };
                EXPECT_EQ(resolver.m_basin_method, fs::mst_method::boruvka);
                EXPECT_EQ(resolver.m_route_method, fs::mst_route_method::basic);
            }
        }

        TEST_P(mst_sink_resolver, resolve_basic)
        {
            auto resolver_ptr
                = std::make_shared<fs::mst_sink_resolver>(GetParam(), fs::mst_route_method::basic);

            auto graph = flow_graph_type(grid, { fs::single_flow_router(), resolver_ptr });

            const auto& new_elevation = graph.update_routes(elevation);

            {
                SCOPED_TRACE("test filled elevation (tiny slope)");

                xt::xtensor<bool, 2> mask_not_filled = xt::ones<bool>(elevation.shape());
                mask_not_filled(1, 3) = false;
                mask_not_filled(2, 2) = false;

                auto new_unchanged = xt::where(mask_not_filled, new_elevation, 0.0);
                auto unchanged = xt::where(mask_not_filled, elevation, 0.0);
                EXPECT_TRUE(xt::all(xt::equal(new_unchanged, unchanged)));

                // pit node (1, 3) is elevated first above the pass (3, 1)
                // then sink node (2, 2) is elevation above the pit node
                EXPECT_TRUE(new_elevation(1, 3) > new_elevation(3, 1));
                EXPECT_TRUE(new_elevation(2, 2) > new_elevation(1, 3));
            }

            {
                SCOPED_TRACE("test resolved flow");

                EXPECT_EQ(graph.impl().receivers()(8, 0), 16);
                EXPECT_EQ(graph.impl().receivers_distance()(8, 0),
                          std::numeric_limits<double>::max());
                EXPECT_EQ(graph.impl().receivers()(12, 0), 8);
                EXPECT_EQ(graph.impl().receivers_distance()(12, 0), std::sqrt(2.0));
            }
        }


        TEST_P(mst_sink_resolver, resolve_carve)
        {
            auto resolver_ptr
                = std::make_shared<fs::mst_sink_resolver>(GetParam(), fs::mst_route_method::carve);

            auto graph = flow_graph_type(grid, { fs::single_flow_router(), resolver_ptr });

            const auto& new_elevation = graph.update_routes(elevation);

            {
                SCOPED_TRACE("test filled elevation (tiny slope)");

                xt::xtensor<bool, 2> mask_not_filled = xt::ones<bool>(elevation.shape());
                mask_not_filled(1, 3) = false;
                mask_not_filled(2, 2) = false;

                auto new_unchanged = xt::where(mask_not_filled, new_elevation, 0.0);
                auto unchanged = xt::where(mask_not_filled, elevation, 0.0);
                EXPECT_TRUE(xt::all(xt::equal(new_unchanged, unchanged)));

                EXPECT_TRUE(new_elevation(2, 2) > new_elevation(3, 1));
                EXPECT_TRUE(new_elevation(1, 3) > new_elevation(2, 2));
            }

            {
                SCOPED_TRACE("test resolved flow");

                EXPECT_EQ(graph.impl().receivers()(8, 0), 12);
                EXPECT_EQ(graph.impl().receivers_distance()(8, 0), std::sqrt(2.0));
                EXPECT_EQ(graph.impl().receivers()(12, 0), 16);
                EXPECT_EQ(graph.impl().receivers_distance()(12, 0), std::sqrt(2.0));
            }
        }

        INSTANTIATE_TEST_SUITE_P(mst_methods,
                                 mst_sink_resolver,
                                 ::testing::Values(fs::mst_method::kruskal,
                                                   fs::mst_method::boruvka));


        /*
         * High-level tests.
         */
        class sink_resolver_extra : public ::testing::Test
        {
        public:
            using grid_type = fs::raster_grid;
            using flow_graph_type = fs::flow_graph<grid_type>;
            using size_type = typename grid_type::size_type;

            grid_type grid = grid_type({ 101, 101 }, { 1.0, 1.0 }, fs::node_status::fixed_value);

            xt::xarray<double> elevation = xt::random::rand<double>(grid.shape());

            /*
             * Check the conservation of drainage area (sum of drainage area at
             * domain fixed boundaries = total domain area).
             */
            template <class SR>
            void test_conservation_drainage_area(const std::shared_ptr<SR>& resolver_ptr)
            {
                using operators_type = fs::flow_operator_sequence<flow_graph_type::impl_type>;
                operators_type ops;

                if constexpr (std::is_same_v<SR, fs::pflood_sink_resolver>)
                {
                    ops = std::move(operators_type(resolver_ptr, fs::single_flow_router()));
                }
                else
                {
                    ops = std::move(operators_type(fs::single_flow_router(), resolver_ptr));
                }

                auto graph = flow_graph_type(grid, std::move(ops));
                graph.update_routes(elevation);
                auto drainage_area = graph.accumulate(1.0);

                double actual = xt::sum(xt::view(drainage_area, xt::keep(0, -1), xt::all()))();
                actual += xt::sum(xt::view(drainage_area, xt::range(1, -1), xt::keep(0, -1)))();

                // assumes grid cell (uniform) area is 1
                double expected = double(elevation.size());

                EXPECT_EQ(actual, expected);
            }

            /*
             * Check that the total number of basins equals the number of (fixed
             * value) boundary nodes at grid borders (i.e., only outer basins).
             */
            template <class SR>
            void test_nb_of_basins(const std::shared_ptr<SR>& resolver_ptr)
            {
                using operators_type = fs::flow_operator_sequence<flow_graph_type::impl_type>;
                operators_type ops;

                if constexpr (std::is_same_v<SR, fs::pflood_sink_resolver>)
                {
                    ops = std::move(operators_type(resolver_ptr, fs::single_flow_router()));
                }
                else
                {
                    ops = std::move(operators_type(fs::single_flow_router(), resolver_ptr));
                }

                auto graph = flow_graph_type(grid, std::move(ops));
                graph.update_routes(elevation);
                auto basins = graph.basins();

                size_type nbasins = xt::amax(basins)() + 1;
                size_type n_border_nodes = (grid.shape()[0] + grid.shape()[1]) * 2 - 4;

                EXPECT_EQ(nbasins, n_border_nodes);
            }
        };


        TEST_F(sink_resolver_extra, pflood_sink_resolver)
        {
            auto resolver_ptr = std::make_shared<pflood_sink_resolver>();

            test_conservation_drainage_area(resolver_ptr);
            test_nb_of_basins(resolver_ptr);
        }


        class mst_sink_resolver_extra
            : public sink_resolver_extra
            , public ::testing::WithParamInterface<fs::mst_sink_resolver>
        {
        };


        TEST_P(mst_sink_resolver_extra, test)
        {
            auto resolver_ptr = std::make_shared<fs::mst_sink_resolver>(GetParam());

            test_conservation_drainage_area(resolver_ptr);
            test_nb_of_basins(resolver_ptr);
        }


        INSTANTIATE_TEST_SUITE_P(
            mst_resolvers,
            mst_sink_resolver_extra,
            ::testing::Values(
                fs::mst_sink_resolver{ fs::mst_method::kruskal, fs::mst_route_method::basic },
                fs::mst_sink_resolver{ fs::mst_method::kruskal, fs::mst_route_method::carve },
                fs::mst_sink_resolver{ fs::mst_method::boruvka, fs::mst_route_method::basic },
                fs::mst_sink_resolver{ fs::mst_method::boruvka, fs::mst_route_method::carve }));
    }
}
