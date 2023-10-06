#include <array>

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

        /*
         * Test single flow routing on a profile grid
         */
        class single_flow_router__profile : public ::testing::Test
        {
        protected:
            using grid_type = fs::profile_grid;
            using size_type = typename grid_type::size_type;
            using flow_graph_type = fs::flow_graph<grid_type>;

            size_type n = static_cast<size_type>(8);

            grid_type grid = grid_type(n, 1., fs::node_status::fixed_value);

            // 3rd node is a pit, 1st and last nodes are base levels
            xt::xtensor<double, 1> elevation{ 0.0, 0.2, 0.1, 0.2, 0.4, 0.6, 0.3, 0.0 };

            flow_graph_type graph = flow_graph_type(grid, { fs::single_flow_router() });

            std::array<size_type, 2> receivers_shape{ grid.size(), 1 };

            void update()
            {
                graph.update_routes(elevation);
            }
        };

        TEST(single_flow_router, ctor)
        {
            EXPECT_EQ(fs::single_flow_router().name(), "single_flow_router");
        }

        TEST_F(single_flow_router__profile, receivers)
        {
            update();

            EXPECT_TRUE(xt::same_shape(graph.impl().receivers().shape(), receivers_shape));

            auto actual = xt::col(graph.impl().receivers(), 0);
            xt::xtensor<size_type, 1> expected{ 0, 0, 2, 2, 3, 6, 7, 7 };

            EXPECT_EQ(actual, expected);
        }

        TEST_F(single_flow_router__profile, receivers_distance)
        {
            update();

            EXPECT_TRUE(xt::same_shape(graph.impl().receivers_distance().shape(), receivers_shape));

            auto actual = xt::col(graph.impl().receivers_distance(), 0);
            xt::xtensor<double, 1> expected{ 0., 1., 0., 1., 1., 1., 1., 0. };

            EXPECT_EQ(actual, expected);
        }

        TEST_F(single_flow_router__profile, receivers_count)
        {
            update();

            auto actual = graph.impl().receivers_count();
            auto expected = xt::ones<std::uint8_t>({ grid.size() });

            EXPECT_EQ(actual, expected);
        }

        TEST_F(single_flow_router__profile, receivers_weight)
        {
            update();

            EXPECT_TRUE(xt::same_shape(graph.impl().receivers_weight().shape(), receivers_shape));

            auto actual = xt::col(graph.impl().receivers_weight(), 0);
            auto expected = xt::ones<double>({ grid.size() });

            EXPECT_EQ(actual, expected);
        }

        TEST_F(single_flow_router__profile, donors)
        {
            update();

            auto actual = graph.impl().donors();
            xt::xarray<int> expected{
                { 1, -1, -1 },  { -1, -1, -1 }, { 2, 3, -1 },  { 4, -1, -1 },
                { -1, -1, -1 }, { -1, -1, -1 }, { 5, -1, -1 }, { 6, -1, -1 }
            };

            EXPECT_EQ(actual, expected);
        }

        TEST_F(single_flow_router__profile, donors_count)
        {
            update();

            auto actual = graph.impl().donors_count();
            xt::xtensor<std::uint8_t, 1> expected{ 1, 0, 2, 1, 0, 0, 1, 1 };

            EXPECT_EQ(actual, expected);
        }

        TEST_F(single_flow_router__profile, dfs_indices)
        {
            update();

            auto actual = graph.impl().dfs_indices();
            xt::xtensor<std::uint8_t, 1> expected{ 0, 1, 2, 3, 4, 7, 6, 5 };

            EXPECT_EQ(actual, expected);
        }

        TEST_F(single_flow_router__profile, mask)
        {
            // mask 4th node only
            xt::xtensor<bool, 1> mask{ false, false, false, true, false, false, false, false };
            graph.set_mask(mask);
            update();

            // masked node has itself as receiver (ignored) and 5th node becomes a pit
            auto actual = xt::col(graph.impl().receivers(), 0);
            xt::xtensor<size_type, 1> expected{ 0, 0, 2, 3, 4, 6, 7, 7 };
            EXPECT_EQ(actual, expected);
        }

        /*
         * Test single flow routing on a raster grid
         */
        template <fs::raster_connect RC>
        class single_flow_router__raster_base : public ::testing::Test
        {
        protected:
            using size_type = typename fs::raster_grid::size_type;

            size_type n = static_cast<size_type>(4);
            std::array<size_type, 2> shape{ { n, n } };

            double dia = std::sqrt(1.1 * 1.1 + 1.2 * 1.2);

            using grid_type = fs::raster_grid_xt<fs::xt_selector, RC>;
            using flow_graph_type = fs::flow_graph<grid_type>;

            // bottom border base-level
            fs::node_status fixed = fs::node_status::fixed_value;
            fs::node_status core = fs::node_status::core;
            fs::node_status looped = fs::node_status::looped;
            fs::raster_boundary_status bottom_base_level{ { core, core, core, fixed } };
            fs::raster_boundary_status bottom_base_level_looped{ { looped, looped, core, fixed } };

            grid_type fixed_grid = grid_type(shape, { 1.1, 1.2 }, bottom_base_level);

            grid_type looped_grid = grid_type(shape, { 1.1, 1.2 }, bottom_base_level_looped);

            flow_graph_type fixed_graph = flow_graph_type(fixed_grid, { fs::single_flow_router() });
            flow_graph_type looped_graph
                = flow_graph_type(looped_grid, { fs::single_flow_router() });

            // planar surface tilted along the y-axis + small carved channel
            // going towards the left border to test looped boundaries
            // (+ small perturbations to avoid ambiguous routing for the bishop case)
            xt::xtensor<double, 2> elevation{ { 0.60, 0.60, 0.61, 0.61 },
                                              { 0.40, 0.40, 0.41, 0.41 },
                                              { 0.05, 0.20, 0.21, 0.21 },
                                              { 0.02, 0.00, 0.11, 0.11 } };

            std::array<size_type, 2> receivers_shape{ fixed_grid.size(), 1 };

            void update()
            {
                fixed_graph.update_routes(elevation);
                looped_graph.update_routes(elevation);
            }
        };


        class single_flow_router__raster_queen
            : public single_flow_router__raster_base<fs::raster_connect::queen>
        {
        protected:
            // receiver / donors all columns except the 1st
            size_type n0plus_cols = 7;
        };


        TEST_F(single_flow_router__raster_queen, receivers)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                EXPECT_TRUE(
                    xt::same_shape(fixed_graph.impl().receivers().shape(), receivers_shape));

                auto actual = xt::col(fixed_graph.impl().receivers(), 0);
                xt::xtensor<size_type, 1> expected{ 4,  5,  6,  7,  8,  8,  10, 11,
                                                    13, 13, 13, 15, 12, 13, 14, 15 };

                EXPECT_EQ(actual, expected);
            }

            {
                SCOPED_TRACE("looped");

                EXPECT_TRUE(
                    xt::same_shape(looped_graph.impl().receivers().shape(), receivers_shape));

                auto actual = xt::col(looped_graph.impl().receivers(), 0);
                xt::xtensor<size_type, 1> expected{ 4,  5,  6,  7, 8,  8,  10, 8,
                                                    13, 13, 13, 8, 12, 13, 14, 15 };

                EXPECT_EQ(actual, expected);
            }
        }

        TEST_F(single_flow_router__raster_queen, receivers_distance)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                EXPECT_TRUE(xt::same_shape(fixed_graph.impl().receivers_distance().shape(),
                                           receivers_shape));

                auto actual = xt::col(fixed_graph.impl().receivers_distance(), 0);
                xt::xtensor<double, 1> expected{ 1.1, 1.1, 1.1, 1.1, 1.1, dia, 1.1, 1.1,
                                                 dia, 1.1, dia, 1.1, 0.0, 0.0, 0.0, 0.0 };

                EXPECT_TRUE(xt::allclose(actual, expected));
            }

            {
                SCOPED_TRACE("looped");

                EXPECT_TRUE(xt::same_shape(looped_graph.impl().receivers_distance().shape(),
                                           receivers_shape));

                auto actual = xt::col(looped_graph.impl().receivers_distance(), 0);
                xt::xtensor<double, 1> expected{ 1.1, 1.1, 1.1, 1.1, 1.1, dia, 1.1, dia,
                                                 dia, 1.1, dia, 1.2, 0.0, 0.0, 0.0, 0.0 };

                EXPECT_TRUE(xt::allclose(actual, expected));
            }
        }

        TEST_F(single_flow_router__raster_queen, receivers_count)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                auto actual = fixed_graph.impl().receivers_count();
                auto expected = xt::ones<std::uint8_t>({ fixed_grid.size() });

                EXPECT_EQ(actual, expected);
            }

            {
                SCOPED_TRACE("looped");

                auto actual = looped_graph.impl().receivers_count();
                auto expected = xt::ones<std::uint8_t>({ looped_grid.size() });

                EXPECT_EQ(actual, expected);
            }
        }

        TEST_F(single_flow_router__raster_queen, receivers_weight)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                EXPECT_TRUE(
                    xt::same_shape(fixed_graph.impl().receivers_weight().shape(), receivers_shape));

                auto actual = xt::col(fixed_graph.impl().receivers_weight(), 0);
                auto expected = xt::ones<double>({ fixed_grid.size() });

                EXPECT_TRUE(xt::allclose(actual, expected));
            }

            {
                SCOPED_TRACE("looped");

                EXPECT_TRUE(xt::same_shape(looped_graph.impl().receivers_weight().shape(),
                                           receivers_shape));

                auto actual = xt::col(looped_graph.impl().receivers_weight(), 0);
                auto expected = xt::ones<double>({ looped_grid.size() });

                EXPECT_TRUE(xt::allclose(actual, expected));
            }
        }

        TEST_F(single_flow_router__raster_queen, donors)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                auto actual = fixed_graph.impl().donors();
                xt::xtensor<int, 2> expected{
                    { -1, -1, -1, -1, -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1, -1, -1, -1, -1 },
                    { -1, -1, -1, -1, -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1, -1, -1, -1, -1 },
                    { 0, -1, -1, -1, -1, -1, -1, -1, -1 },  { 1, -1, -1, -1, -1, -1, -1, -1, -1 },
                    { 2, -1, -1, -1, -1, -1, -1, -1, -1 },  { 3, -1, -1, -1, -1, -1, -1, -1, -1 },
                    { 4, 5, -1, -1, -1, -1, -1, -1, -1 },   { -1, -1, -1, -1, -1, -1, -1, -1, -1 },
                    { 6, -1, -1, -1, -1, -1, -1, -1, -1 },  { 7, -1, -1, -1, -1, -1, -1, -1, -1 },
                    { -1, -1, -1, -1, -1, -1, -1, -1, -1 }, { 8, 9, 10, -1, -1, -1, -1, -1, -1 },
                    { -1, -1, -1, -1, -1, -1, -1, -1, -1 }, { 11, -1, -1, -1, -1, -1, -1, -1, -1 }
                };

                EXPECT_EQ(actual, expected);
            }

            {
                SCOPED_TRACE("looped");

                auto actual = looped_graph.impl().donors();
                xt::xtensor<int, 2> expected{
                    { -1, -1, -1, -1, -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1, -1, -1, -1, -1 },
                    { -1, -1, -1, -1, -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1, -1, -1, -1, -1 },
                    { 0, -1, -1, -1, -1, -1, -1, -1, -1 },  { 1, -1, -1, -1, -1, -1, -1, -1, -1 },
                    { 2, -1, -1, -1, -1, -1, -1, -1, -1 },  { 3, -1, -1, -1, -1, -1, -1, -1, -1 },
                    { 4, 5, 7, 11, -1, -1, -1, -1, -1 },    { -1, -1, -1, -1, -1, -1, -1, -1, -1 },
                    { 6, -1, -1, -1, -1, -1, -1, -1, -1 },  { -1, -1, -1, -1, -1, -1, -1, -1, -1 },
                    { -1, -1, -1, -1, -1, -1, -1, -1, -1 }, { 8, 9, 10, -1, -1, -1, -1, -1, -1 },
                    { -1, -1, -1, -1, -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1, -1, -1, -1, -1 }
                };

                EXPECT_EQ(actual, expected);
            }
        }

        TEST_F(single_flow_router__raster_queen, donors_count)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                auto actual = fixed_graph.impl().donors_count();
                xt::xtensor<std::uint8_t, 1> expected{ 0, 0, 0, 0, 1, 1, 1, 1,
                                                       2, 0, 1, 1, 0, 3, 0, 1 };

                EXPECT_EQ(actual, expected);
            }

            {
                SCOPED_TRACE("looped");

                auto actual = looped_graph.impl().donors_count();
                xt::xtensor<std::uint8_t, 1> expected{ 0, 0, 0, 0, 1, 1, 1, 1,
                                                       4, 0, 1, 0, 0, 3, 0, 0 };

                EXPECT_EQ(actual, expected);
            }
        }

        TEST_F(single_flow_router__raster_queen, dfs_indices)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                auto actual = fixed_graph.impl().dfs_indices();
                xt::xtensor<std::uint8_t, 1> expected{ 12, 13, 8, 9,  10, 6,  2, 4,
                                                       5,  1,  0, 14, 15, 11, 7, 3 };

                EXPECT_EQ(actual, expected);
            }

            {
                SCOPED_TRACE("looped");

                auto actual = looped_graph.impl().dfs_indices();
                xt::xtensor<std::uint8_t, 1> expected{ 12, 13, 8,  9, 10, 6, 2,  4,
                                                       5,  7,  11, 3, 1,  0, 14, 15 };

                EXPECT_EQ(actual, expected);
            }
        }


        class single_flow_router__raster_rook
            : public single_flow_router__raster_base<fs::raster_connect::rook>
        {
        protected:
            // receiver / donors all columns except the 1st
            size_type n0plus_cols = 3;
        };


        TEST_F(single_flow_router__raster_rook, receivers)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                EXPECT_TRUE(
                    xt::same_shape(fixed_graph.impl().receivers().shape(), receivers_shape));

                auto actual = xt::col(fixed_graph.impl().receivers(), 0);
                xt::xtensor<size_type, 1> expected{ 4,  5,  6,  7,  8,  9,  10, 11,
                                                    12, 13, 14, 15, 12, 13, 14, 15 };

                EXPECT_EQ(actual, expected);
            }

            {
                SCOPED_TRACE("looped");

                EXPECT_TRUE(
                    xt::same_shape(looped_graph.impl().receivers().shape(), receivers_shape));

                auto actual = xt::col(looped_graph.impl().receivers(), 0);
                xt::xtensor<size_type, 1> expected{ 4,  5,  6,  7, 8,  9,  10, 11,
                                                    12, 13, 14, 8, 12, 13, 14, 15 };

                EXPECT_EQ(actual, expected);
            }
        }

        TEST_F(single_flow_router__raster_rook, receivers_distance)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                EXPECT_TRUE(xt::same_shape(fixed_graph.impl().receivers_distance().shape(),
                                           receivers_shape));

                auto actual = xt::col(fixed_graph.impl().receivers_distance(), 0);
                xt::xtensor<double, 1> expected{ 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1,
                                                 1.1, 1.1, 1.1, 1.1, 0.0, 0.0, 0.0, 0.0 };

                EXPECT_TRUE(xt::allclose(actual, expected));
            }
            {
                SCOPED_TRACE("looped");

                EXPECT_TRUE(xt::same_shape(looped_graph.impl().receivers_distance().shape(),
                                           receivers_shape));

                auto actual = xt::col(looped_graph.impl().receivers_distance(), 0);
                xt::xtensor<double, 1> expected{ 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1,
                                                 1.1, 1.1, 1.1, 1.2, 0.0, 0.0, 0.0, 0.0 };

                EXPECT_TRUE(xt::allclose(actual, expected));
            }
        }

        TEST_F(single_flow_router__raster_rook, receivers_count)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                auto actual = fixed_graph.impl().receivers_count();
                auto expected = xt::ones<std::uint8_t>({ fixed_grid.size() });

                EXPECT_EQ(actual, expected);
            }

            {
                SCOPED_TRACE("looped");

                auto actual = looped_graph.impl().receivers_count();
                auto expected = xt::ones<std::uint8_t>({ looped_grid.size() });

                EXPECT_EQ(actual, expected);
            }
        }

        TEST_F(single_flow_router__raster_rook, receivers_weight)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                EXPECT_TRUE(
                    xt::same_shape(fixed_graph.impl().receivers_weight().shape(), receivers_shape));

                auto actual = xt::col(fixed_graph.impl().receivers_weight(), 0);
                auto expected = xt::ones<double>({ fixed_grid.size() });

                EXPECT_TRUE(xt::allclose(actual, expected));
            }

            {
                SCOPED_TRACE("looped");

                EXPECT_TRUE(xt::same_shape(looped_graph.impl().receivers_weight().shape(),
                                           receivers_shape));

                auto actual = xt::col(looped_graph.impl().receivers_weight(), 0);
                auto expected = xt::ones<double>({ looped_grid.size() });

                EXPECT_TRUE(xt::allclose(actual, expected));
            }
        }

        TEST_F(single_flow_router__raster_rook, donors)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                auto actual = fixed_graph.impl().donors();
                xt::xtensor<int, 2> expected{ { -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1 },
                                              { -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1 },
                                              { 0, -1, -1, -1, -1 },  { 1, -1, -1, -1, -1 },
                                              { 2, -1, -1, -1, -1 },  { 3, -1, -1, -1, -1 },
                                              { 4, -1, -1, -1, -1 },  { 5, -1, -1, -1, -1 },
                                              { 6, -1, -1, -1, -1 },  { 7, -1, -1, -1, -1 },
                                              { 8, -1, -1, -1, -1 },  { 9, -1, -1, -1, -1 },
                                              { 10, -1, -1, -1, -1 }, { 11, -1, -1, -1, -1 } };

                EXPECT_EQ(actual, expected);
            }

            {
                SCOPED_TRACE("looped");

                auto actual = looped_graph.impl().donors();
                xt::xtensor<int, 2> expected{ { -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1 },
                                              { -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1 },
                                              { 0, -1, -1, -1, -1 },  { 1, -1, -1, -1, -1 },
                                              { 2, -1, -1, -1, -1 },  { 3, -1, -1, -1, -1 },
                                              { 4, 11, -1, -1, -1 },  { 5, -1, -1, -1, -1 },
                                              { 6, -1, -1, -1, -1 },  { 7, -1, -1, -1, -1 },
                                              { 8, -1, -1, -1, -1 },  { 9, -1, -1, -1, -1 },
                                              { 10, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1 } };

                EXPECT_EQ(actual, expected);
            }
        }

        TEST_F(single_flow_router__raster_rook, donors_count)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                auto actual = fixed_graph.impl().donors_count();
                xt::xtensor<std::uint8_t, 1> expected{ 0, 0, 0, 0, 1, 1, 1, 1,
                                                       1, 1, 1, 1, 1, 1, 1, 1 };

                EXPECT_EQ(actual, expected);
            }

            {
                SCOPED_TRACE("looped");

                auto actual = looped_graph.impl().donors_count();
                xt::xtensor<std::uint8_t, 1> expected{ 0, 0, 0, 0, 1, 1, 1, 1,
                                                       2, 1, 1, 1, 1, 1, 1, 0 };

                EXPECT_EQ(actual, expected);
            }
        }

        TEST_F(single_flow_router__raster_rook, dfs_indices)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                auto actual = fixed_graph.impl().dfs_indices();
                xt::xtensor<std::uint8_t, 1> expected{ 12, 8,  4, 0, 13, 9,  5, 1,
                                                       14, 10, 6, 2, 15, 11, 7, 3 };

                EXPECT_EQ(actual, expected);
            }

            {
                SCOPED_TRACE("looped");

                auto actual = looped_graph.impl().dfs_indices();
                xt::xtensor<std::uint8_t, 1> expected{ 12, 8, 4, 11, 7,  3, 0, 13,
                                                       9,  5, 1, 14, 10, 6, 2, 15 };

                EXPECT_EQ(actual, expected);
            }
        }


        class single_flow_router__raster_bishop
            : public single_flow_router__raster_base<fs::raster_connect::bishop>
        {
        protected:
            // receiver / donors all columns except the 1st
            size_type n0plus_cols = 3;
        };


        TEST_F(single_flow_router__raster_bishop, receivers)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                EXPECT_TRUE(
                    xt::same_shape(fixed_graph.impl().receivers().shape(), receivers_shape));

                auto actual = xt::col(fixed_graph.impl().receivers(), 0);
                xt::xtensor<size_type, 1> expected{ 5,  4,  5,  6,  9,  8,  9,  10,
                                                    13, 12, 13, 14, 12, 13, 14, 15 };

                EXPECT_EQ(actual, expected);
            }

            {
                SCOPED_TRACE("looped");

                EXPECT_TRUE(
                    xt::same_shape(looped_graph.impl().receivers().shape(), receivers_shape));

                auto actual = xt::col(looped_graph.impl().receivers(), 0);
                xt::xtensor<size_type, 1> expected{ 5,  4,  5,  4,  9,  8,  9,  8,
                                                    13, 12, 13, 12, 12, 13, 14, 15 };

                EXPECT_EQ(actual, expected);
            }
        }

        TEST_F(single_flow_router__raster_bishop, receivers_distance)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                EXPECT_TRUE(xt::same_shape(fixed_graph.impl().receivers_distance().shape(),
                                           receivers_shape));

                auto actual = xt::col(fixed_graph.impl().receivers_distance(), 0);
                xt::xtensor<double, 1> expected{ dia, dia, dia, dia, dia, dia, dia, dia,
                                                 dia, dia, dia, dia, 0.0, 0.0, 0.0, 0.0 };

                EXPECT_TRUE(xt::allclose(actual, expected));
            }

            {
                SCOPED_TRACE("looped");

                EXPECT_TRUE(xt::same_shape(looped_graph.impl().receivers_distance().shape(),
                                           receivers_shape));

                auto actual = xt::col(looped_graph.impl().receivers_distance(), 0);
                xt::xtensor<double, 1> expected{ dia, dia, dia, dia, dia, dia, dia, dia,
                                                 dia, dia, dia, dia, 0.0, 0.0, 0.0, 0.0 };

                EXPECT_TRUE(xt::allclose(actual, expected));
            }
        }

        TEST_F(single_flow_router__raster_bishop, receivers_count)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                auto actual = fixed_graph.impl().receivers_count();
                auto expected = xt::ones<std::uint8_t>({ fixed_grid.size() });

                EXPECT_EQ(actual, expected);
            }

            {
                SCOPED_TRACE("looped");

                auto actual = looped_graph.impl().receivers_count();
                auto expected = xt::ones<std::uint8_t>({ looped_grid.size() });

                EXPECT_EQ(actual, expected);
            }
        }

        TEST_F(single_flow_router__raster_bishop, receivers_weight)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                EXPECT_TRUE(
                    xt::same_shape(fixed_graph.impl().receivers_weight().shape(), receivers_shape));

                auto actual = xt::col(fixed_graph.impl().receivers_weight(), 0);
                auto expected = xt::ones<double>({ fixed_grid.size() });

                EXPECT_TRUE(xt::allclose(actual, expected));
            }

            {
                SCOPED_TRACE("looped");

                EXPECT_TRUE(xt::same_shape(looped_graph.impl().receivers_weight().shape(),
                                           receivers_shape));

                auto actual = xt::col(looped_graph.impl().receivers_weight(), 0);
                auto expected = xt::ones<double>({ looped_grid.size() });

                EXPECT_TRUE(xt::allclose(actual, expected));
            }
        }

        TEST_F(single_flow_router__raster_bishop, donors)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                auto actual = fixed_graph.impl().donors();
                xt::xtensor<int, 2> expected{ { -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1 },
                                              { -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1 },
                                              { 1, -1, -1, -1, -1 },  { 0, 2, -1, -1, -1 },
                                              { 3, -1, -1, -1, -1 },  { -1, -1, -1, -1, -1 },
                                              { 5, -1, -1, -1, -1 },  { 4, 6, -1, -1, -1 },
                                              { 7, -1, -1, -1, -1 },  { -1, -1, -1, -1, -1 },
                                              { 9, -1, -1, -1, -1 },  { 8, 10, -1, -1, -1 },
                                              { 11, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1 } };

                EXPECT_EQ(actual, expected);
            }

            {
                SCOPED_TRACE("looped");

                auto actual = looped_graph.impl().donors();
                xt::xtensor<int, 2> expected{ { -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1 },
                                              { -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1 },
                                              { 1, 3, -1, -1, -1 },   { 0, 2, -1, -1, -1 },
                                              { -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1 },
                                              { 5, 7, -1, -1, -1 },   { 4, 6, -1, -1, -1 },
                                              { -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1 },
                                              { 9, 11, -1, -1, -1 },  { 8, 10, -1, -1, -1 },
                                              { -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1 } };

                EXPECT_EQ(actual, expected);
            }
        }

        TEST_F(single_flow_router__raster_bishop, donors_count)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                auto actual = fixed_graph.impl().donors_count();
                xt::xtensor<std::uint8_t, 1> expected{ 0, 0, 0, 0, 1, 2, 1, 0,
                                                       1, 2, 1, 0, 1, 2, 1, 0 };

                EXPECT_EQ(actual, expected);
            }

            {
                SCOPED_TRACE("looped");

                auto actual = looped_graph.impl().donors_count();
                xt::xtensor<std::uint8_t, 1> expected{ 0, 0, 0, 0, 2, 2, 0, 0,
                                                       2, 2, 0, 0, 2, 2, 0, 0 };

                EXPECT_EQ(actual, expected);
            }
        }

        TEST_F(single_flow_router__raster_bishop, dfs_indices)
        {
            update();

            {
                SCOPED_TRACE("non-looped");

                auto actual = fixed_graph.impl().dfs_indices();
                xt::xtensor<std::uint8_t, 1> expected{ 12, 9, 4, 6, 3, 1,  13, 8,
                                                       10, 7, 5, 0, 2, 14, 11, 15 };

                EXPECT_EQ(actual, expected);
            }

            {
                SCOPED_TRACE("looped");

                auto actual = looped_graph.impl().dfs_indices();
                xt::xtensor<std::uint8_t, 1> expected{ 12, 9,  11, 4, 6, 1, 3,  13,
                                                       8,  10, 5,  7, 0, 2, 14, 15 };

                EXPECT_EQ(actual, expected);
            }
        }
    }
}
