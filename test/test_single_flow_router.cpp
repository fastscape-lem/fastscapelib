#include "fastscapelib/flow_graph.hpp"
#include "fastscapelib/flow_router.hpp"
#include "fastscapelib/flow_router_factory.hpp"
#include "fastscapelib/sink_resolver.hpp"
#include "fastscapelib/grid/raster_grid.hpp"

#include "xtensor/xtensor.hpp"
#include "xtensor/xio.hpp"

#include "gtest/gtest.h"


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace testing
    {

        class single_flow_router__profile : public ::testing::Test
        {
        protected:
            using grid_type = fs::profile_grid;
            using size_type = typename grid_type::size_type;
            using flow_graph_type = fs::flow_graph<grid_type, double>;

            size_type n = static_cast<size_type>(8);

            grid_type fixed_grid = grid_type(n, 1., fs::node_status::fixed_value_boundary);
            grid_type looped_grid = grid_type(n, 1., fs::node_status::looped_boundary);

            xt::xtensor<double, 1> elevation{ 0.82, 0.16, 0.14, 0.20, 0.71, 0.97, 0.41, 0.09 };

            flow_graph_type fixed_graph
                = flow_graph_type(fixed_grid,
                                  std::make_unique<fs::single_flow_router<flow_graph_type>>(),
                                  std::make_unique<fs::no_sink_resolver<flow_graph_type>>());

            flow_graph_type looped_graph
                = flow_graph_type(looped_grid,
                                  std::make_unique<fs::single_flow_router<flow_graph_type>>(),
                                  std::make_unique<fs::no_sink_resolver<flow_graph_type>>());

            void update()
            {
                fixed_graph.update_routes(elevation);
                looped_graph.update_routes(elevation);
            }
        };

        TEST_F(single_flow_router__profile, receivers)
        {
            update();

            xt::xtensor<size_type, 1> expected_fixed_receivers{ 1, 2, 2, 2, 3, 6, 7, 7 };
            EXPECT_TRUE(
                xt::all(xt::equal(xt::col(fixed_graph.receivers(), 0), expected_fixed_receivers)));
            EXPECT_TRUE(xt::all(xt::equal(xt::col(fixed_graph.receivers(), 1),
                                          xt::ones_like(expected_fixed_receivers) * -1)));

            xt::xtensor<size_type, 1> expected_looped_receivers{ 7, 2, 2, 2, 3, 6, 7, 7 };
            EXPECT_TRUE(xt::all(
                xt::equal(xt::col(looped_graph.receivers(), 0), expected_looped_receivers)));
            EXPECT_TRUE(xt::all(xt::equal(xt::col(looped_graph.receivers(), 1),
                                          xt::ones_like(expected_fixed_receivers) * -1)));
        }

        TEST_F(single_flow_router__profile, receivers_distance)
        {
            update();

            xt::xtensor<double, 1> expected_fixed_receivers_distance{
                1., 1., 0., 1., 1., 1., 1., 0.
            };
            EXPECT_TRUE(xt::all(xt::equal(xt::col(fixed_graph.receivers_distance(), 0),
                                          expected_fixed_receivers_distance)));
            EXPECT_TRUE(xt::all(xt::equal(xt::col(fixed_graph.receivers_distance(), 1),
                                          xt::ones_like(expected_fixed_receivers_distance) * -1)));

            xt::xtensor<double, 1> expected_looped_receivers_distance{ 1., 1., 0., 1.,
                                                                       1., 1., 1., 0. };
            EXPECT_TRUE(xt::all(xt::equal(xt::col(looped_graph.receivers_distance(), 0),
                                          expected_looped_receivers_distance)));
            EXPECT_TRUE(xt::all(xt::equal(xt::col(looped_graph.receivers_distance(), 1),
                                          xt::ones_like(expected_fixed_receivers_distance) * -1)));
        }

        TEST_F(single_flow_router__profile, receivers_count)
        {
            update();

            EXPECT_TRUE(
                xt::all(xt::equal(fixed_graph.receivers_count(), xt::ones<std::uint8_t>({ 8 }))));

            EXPECT_TRUE(
                xt::all(xt::equal(looped_graph.receivers_count(), xt::ones<std::uint8_t>({ 8 }))));
        }

        TEST_F(single_flow_router__profile, receivers_weight)
        {
            update();

            EXPECT_TRUE(xt::all(
                xt::equal(xt::col(fixed_graph.receivers_weight(), 0), xt::ones<double>({ 8 }))));
            EXPECT_TRUE(xt::all(
                xt::equal(xt::col(fixed_graph.receivers_weight(), 1), xt::zeros<double>({ 8 }))));

            EXPECT_TRUE(xt::all(
                xt::equal(xt::col(looped_graph.receivers_weight(), 0), xt::ones<double>({ 8 }))));
            EXPECT_TRUE(xt::all(
                xt::equal(xt::col(looped_graph.receivers_weight(), 1), xt::zeros<double>({ 8 }))));
        }

        TEST_F(single_flow_router__profile, donors)
        {
            update();

            xt::xarray<size_type> expected_fixed_donors = xt::ones<size_type>({ 8, 3 }) * -1;
            expected_fixed_donors(1, 0) = 0;
            xt::view(expected_fixed_donors, 2, xt::all()) = xt::xarray<size_type>({ 1, 2, 3 });
            expected_fixed_donors(3, 0) = 4;
            expected_fixed_donors(6, 0) = 5;
            xt::view(expected_fixed_donors, 7, xt::range(0, 2)) = xt::xarray<size_type>({ 6, 7 });

            EXPECT_TRUE(xt::all(xt::equal(fixed_graph.donors(), expected_fixed_donors)));

            xt::xarray<size_type> expected_looped_donors = xt::ones<size_type>({ 8, 3 }) * -1;
            xt::view(expected_looped_donors, 2, xt::all()) = xt::xarray<size_type>({ 1, 2, 3 });
            expected_looped_donors(3, 0) = 4;
            expected_looped_donors(6, 0) = 5;
            xt::view(expected_looped_donors, 7, xt::all()) = xt::xarray<size_type>({ 0, 6, 7 });

            EXPECT_TRUE(xt::all(xt::equal(looped_graph.donors(), expected_looped_donors)));
        }

        TEST_F(single_flow_router__profile, donors_count)
        {
            update();

            xt::xtensor<std::uint8_t, 1> expected_fixed_donors_count{ 0, 1, 3, 1, 0, 0, 1, 2 };
            EXPECT_TRUE(
                xt::all(xt::equal(fixed_graph.donors_count(), expected_fixed_donors_count)));

            xt::xtensor<std::uint8_t, 1> expected_looped_donors_count{ 0, 0, 3, 1, 0, 0, 1, 3 };
            EXPECT_TRUE(
                xt::all(xt::equal(looped_graph.donors_count(), expected_looped_donors_count)));
        }

        TEST_F(single_flow_router__profile, dfs_stack)
        {
            update();

            xt::xtensor<std::uint8_t, 1> expected_fixed_stack{ 2, 1, 0, 3, 4, 7, 6, 5 };
            EXPECT_TRUE(xt::all(xt::equal(fixed_graph.dfs_stack(), expected_fixed_stack)));

            xt::xtensor<std::uint8_t, 1> expected_looped_stack{ 2, 1, 3, 4, 7, 0, 6, 5 };
            EXPECT_TRUE(xt::all(xt::equal(looped_graph.dfs_stack(), expected_looped_stack)));
        }

        TEST_F(single_flow_router__profile, dfs_iterators)
        {
            update();

            xt::xtensor<std::uint8_t, 1> expected_fixed_stack{ 2, 1, 0, 3, 4, 7, 6, 5 };

            auto bfs = fixed_graph.dfs_cbegin();
            for (std::size_t i = 0; i < fixed_graph.size(); ++i)
            {
                EXPECT_EQ(*bfs, expected_fixed_stack[i]);
                ++bfs;
            }

            auto rbfs = fixed_graph.dfs_crbegin();
            for (std::size_t i = fixed_graph.size(); i != 0; --i)
            {
                EXPECT_EQ(*rbfs, expected_fixed_stack[i - 1]);
                ++rbfs;
            }

            std::size_t sum = 0;
            for (auto it = fixed_graph.dfs_cbegin(); it != fixed_graph.dfs_cend(); ++it)
            {
                sum += *it;
            }
            EXPECT_EQ(sum, 28u);

            sum = 0;
            for (auto it = fixed_graph.dfs_crbegin(); it != fixed_graph.dfs_crend(); ++it)
            {
                sum += *it;
            }
            EXPECT_EQ(sum, 28u);
        }


        class single_flow_router__raster_base : public ::testing::Test
        {
        protected:
            using size_type = typename fs::raster_grid::size_type;

            size_type n = static_cast<size_type>(4);
            std::array<size_type, 2> shape{ { n, n } };

            double dia = std::sqrt(1.1 * 1.1 + 1.2 * 1.2);
        };

        class single_flow_router__raster_queen : public single_flow_router__raster_base
        {
        protected:
            using grid_type = fs::raster_grid_xt<fs::xtensor_selector, fs::raster_connect::queen>;
            using flow_graph_type = fs::flow_graph<grid_type, double>;

            grid_type fixed_grid
                = grid_type(shape, { 1.1, 1.2 }, fs::node_status::fixed_value_boundary);
            grid_type looped_grid
                = grid_type(shape, { 1.1, 1.2 }, fs::node_status::looped_boundary);

            flow_graph_type fixed_graph
                = flow_graph_type(fixed_grid,
                                  std::make_unique<fs::single_flow_router<flow_graph_type>>(),
                                  std::make_unique<fs::no_sink_resolver<flow_graph_type>>());

            flow_graph_type looped_graph
                = flow_graph_type(looped_grid,
                                  std::make_unique<fs::single_flow_router<flow_graph_type>>(),
                                  std::make_unique<fs::no_sink_resolver<flow_graph_type>>());

            xt::xtensor<double, 2> elevation{ { 0.82, 0.16, 0.14, 0.20 },
                                              { 0.71, 0.97, 0.41, 0.09 },
                                              { 0.49, 0.01, 0.19, 0.38 },
                                              { 0.29, 0.82, 0.09, 0.88 } };

            void update()
            {
                fixed_graph.update_routes(elevation);
                looped_graph.update_routes(elevation);
            }
        };

        TEST_F(single_flow_router__raster_queen, receivers)
        {
            update();

            xt::xtensor<size_type, 1> expected_fixed_receivers{ 1, 2, 7, 7, 9, 9, 7, 7,
                                                                9, 9, 9, 7, 9, 9, 9, 14 };

            EXPECT_TRUE(
                xt::all(xt::equal(xt::col(fixed_graph.receivers(), 0), expected_fixed_receivers)));
            EXPECT_TRUE(xt::all(xt::equal(xt::view(fixed_graph.receivers(), 0, xt::range(1, 8)),
                                          xt::ones<size_type>({ 16, 7 }) * -1)));

            xt::xtensor<size_type, 1> expected_looped_receivers{ 1, 14, 14, 7, 7, 9, 7, 7,
                                                                 9, 9,  9,  7, 9, 9, 9, 14 };

            EXPECT_TRUE(xt::all(
                xt::equal(xt::col(looped_graph.receivers(), 0), expected_looped_receivers)));
            EXPECT_TRUE(xt::all(xt::equal(xt::view(looped_graph.receivers(), 0, xt::range(1, 8)),
                                          xt::ones<size_type>({ 16, 7 }) * -1)));
        }

        TEST_F(single_flow_router__raster_queen, receivers_distance)
        {
            update();

            xt::xtensor<double, 1> expected_fixed_receivers_distance{ 1.2, 1.2, dia, 1.1, dia, 1.1,
                                                                      1.2, 0.0, 1.2, 0.0, 1.2, 1.1,
                                                                      dia, 1.1, dia, 1.2 };

            EXPECT_TRUE(xt::allclose(xt::col(fixed_graph.receivers_distance(), 0),
                                     expected_fixed_receivers_distance));
            EXPECT_TRUE(
                xt::all(xt::equal(xt::view(fixed_graph.receivers_distance(), 0, xt::range(1, 8)),
                                  xt::ones<double>({ 16, 7 }) * -1)));

            xt::xtensor<double, 1> expected_looped_receivers_distance{ 1.2, dia, 1.1, 1.1, 1.2, 1.1,
                                                                       1.2, 0.0, 1.2, 0.0, 1.2, 1.1,
                                                                       dia, 1.1, dia, 1.2 };

            EXPECT_TRUE(xt::allclose(xt::col(looped_graph.receivers_distance(), 0),
                                     expected_looped_receivers_distance));
            EXPECT_TRUE(
                xt::all(xt::equal(xt::view(looped_graph.receivers_distance(), 0, xt::range(1, 8)),
                                  xt::ones<double>({ 16, 7 }) * -1)));
        }

        TEST_F(single_flow_router__raster_queen, receivers_count)
        {
            update();

            EXPECT_TRUE(
                xt::all(xt::equal(fixed_graph.receivers_count(), xt::ones<std::uint8_t>({ 16 }))));

            EXPECT_TRUE(
                xt::all(xt::equal(looped_graph.receivers_count(), xt::ones<std::uint8_t>({ 16 }))));
        }

        TEST_F(single_flow_router__raster_queen, receivers_weight)
        {
            update();

            EXPECT_TRUE(xt::all(
                xt::equal(xt::col(fixed_graph.receivers_weight(), 0), xt::ones<double>({ 16 }))));
            EXPECT_TRUE(xt::all(
                xt::equal(xt::view(fixed_graph.receivers_weight(), xt::all(), xt::range(1, 8)),
                          xt::zeros<double>({ 16, 7 }))));

            EXPECT_TRUE(xt::all(
                xt::equal(xt::col(looped_graph.receivers_weight(), 0), xt::ones<double>({ 16 }))));
            EXPECT_TRUE(xt::all(
                xt::equal(xt::view(looped_graph.receivers_weight(), xt::all(), xt::range(1, 8)),
                          xt::zeros<double>({ 16, 7 }))));
        }

        TEST_F(single_flow_router__raster_queen, donors)
        {
            update();

            xt::xarray<size_type> expected_fixed_donors = xt::ones<size_type>({ 16, 9 }) * -1;
            expected_fixed_donors(1, 0) = 0;
            expected_fixed_donors(2, 0) = 1;
            xt::view(expected_fixed_donors, 7, xt::range(0, 5))
                = xt::xarray<size_type>({ 2, 3, 6, 7, 11 });
            xt::view(expected_fixed_donors, 9, xt::range(0, 8))
                = xt::xarray<size_type>({ 4, 5, 8, 9, 10, 12, 13, 14 });
            expected_fixed_donors(14, 0) = 15;

            EXPECT_TRUE(xt::all(xt::equal(fixed_graph.donors(), expected_fixed_donors)));

            xt::xarray<size_type> expected_looped_donors = xt::ones<size_type>({ 16, 9 }) * -1;
            expected_looped_donors(1, 0) = 0;
            xt::view(expected_looped_donors, 7, xt::range(0, 5))
                = xt::xarray<size_type>({ 3, 4, 6, 7, 11 });
            xt::view(expected_looped_donors, 9, xt::range(0, 7))
                = xt::xarray<size_type>({ 5, 8, 9, 10, 12, 13, 14 });
            xt::view(expected_looped_donors, 14, xt::range(0, 3))
                = xt::xarray<size_type>({ 1, 2, 15 });

            EXPECT_TRUE(xt::all(xt::equal(looped_graph.donors(), expected_looped_donors)));
        }

        TEST_F(single_flow_router__raster_queen, donors_count)
        {
            update();

            xt::xtensor<std::uint8_t, 1> expected_fixed_donors_count{ 0, 1, 1, 0, 0, 0, 0, 5,
                                                                      0, 8, 0, 0, 0, 0, 1, 0 };

            EXPECT_TRUE(
                xt::all(xt::equal(fixed_graph.donors_count(), expected_fixed_donors_count)));

            xt::xtensor<std::uint8_t, 1> expected_looped_donors_count{ 0, 1, 0, 0, 0, 0, 0, 5,
                                                                       0, 7, 0, 0, 0, 0, 3, 0 };

            EXPECT_TRUE(
                xt::all(xt::equal(looped_graph.donors_count(), expected_looped_donors_count)));
        }


        TEST_F(single_flow_router__raster_queen, dfs_stack)
        {
            update();

            xt::xtensor<std::uint8_t, 1> expected_fixed_stack{ 7, 2, 1, 0,  3,  6,  11, 9,
                                                               4, 5, 8, 10, 12, 13, 14, 15 };
            EXPECT_TRUE(xt::all(xt::equal(fixed_graph.dfs_stack(), expected_fixed_stack)));

            xt::xtensor<std::uint8_t, 1> expected_looped_stack{ 7,  3,  4,  6,  11, 9, 5, 8,
                                                                10, 12, 13, 14, 1,  0, 2, 15 };
            EXPECT_TRUE(xt::all(xt::equal(looped_graph.dfs_stack(), expected_looped_stack)));
        }

        TEST_F(single_flow_router__raster_queen, dfs_iterators)
        {
            update();

            xt::xtensor<std::uint8_t, 1> expected_fixed_stack{ 7, 2, 1, 0,  3,  6,  11, 9,
                                                               4, 5, 8, 10, 12, 13, 14, 15 };

            auto bfs = fixed_graph.dfs_cbegin();
            for (std::size_t i = 0; i < fixed_graph.size(); ++i)
            {
                EXPECT_EQ(*bfs, expected_fixed_stack[i]);
                ++bfs;
            }

            auto rbfs = fixed_graph.dfs_crbegin();
            for (std::size_t i = fixed_graph.size(); i != 0; --i)
            {
                EXPECT_EQ(*rbfs, expected_fixed_stack[i - 1]);
                ++rbfs;
            }

            std::size_t sum = 0;
            for (auto it = fixed_graph.dfs_cbegin(); it != fixed_graph.dfs_cend(); ++it)
            {
                sum += *it;
            }
            EXPECT_EQ(sum, 120u);

            sum = 0;
            for (auto it = fixed_graph.dfs_crbegin(); it != fixed_graph.dfs_crend(); ++it)
            {
                sum += *it;
            }
            EXPECT_EQ(sum, 120u);
        }

        class single_flow_router__raster_rook : public single_flow_router__raster_base
        {
        protected:
            using grid_type = fs::raster_grid_xt<fs::xtensor_selector, fs::raster_connect::rook>;
            using flow_graph_type = fs::flow_graph<grid_type, double>;

            grid_type fixed_grid
                = grid_type(shape, { 1.1, 1.2 }, fs::node_status::fixed_value_boundary);
            grid_type looped_grid
                = grid_type(shape, { 1.1, 1.2 }, fs::node_status::looped_boundary);

            flow_graph_type fixed_graph
                = flow_graph_type(fixed_grid,
                                  std::make_unique<fs::single_flow_router<flow_graph_type>>(),
                                  std::make_unique<fs::no_sink_resolver<flow_graph_type>>());

            flow_graph_type looped_graph
                = flow_graph_type(looped_grid,
                                  std::make_unique<fs::single_flow_router<flow_graph_type>>(),
                                  std::make_unique<fs::no_sink_resolver<flow_graph_type>>());

            xt::xtensor<double, 2> elevation{ { 0.82, 0.16, 0.14, 0.20 },
                                              { 0.71, 0.97, 0.41, 0.09 },
                                              { 0.49, 0.01, 0.10, 0.38 },
                                              { 0.29, 0.82, 0.09, 0.88 } };

            void update()
            {
                fixed_graph.update_routes(elevation);
                looped_graph.update_routes(elevation);
            }
        };

        TEST_F(single_flow_router__raster_rook, receivers)
        {
            update();

            xt::xtensor<size_type, 1> expected_fixed_receivers{ 1, 2, 2, 7, 8,  9, 10, 7,
                                                                9, 9, 9, 7, 12, 9, 14, 14 };

            EXPECT_TRUE(
                xt::all(xt::equal(xt::col(fixed_graph.receivers(), 0), expected_fixed_receivers)));
            EXPECT_TRUE(xt::all(xt::equal(xt::view(fixed_graph.receivers(), 0, xt::range(1, 4)),
                                          xt::ones<size_type>({ 16, 3 }) * -1)));

            xt::xtensor<size_type, 1> expected_looped_receivers{ 1, 2, 14, 7, 7,  9, 10, 7,
                                                                 9, 9, 9,  7, 12, 9, 14, 14 };

            EXPECT_TRUE(xt::all(
                xt::equal(xt::col(looped_graph.receivers(), 0), expected_looped_receivers)));
            EXPECT_TRUE(xt::all(xt::equal(xt::view(looped_graph.receivers(), 0, xt::range(1, 4)),
                                          xt::ones<size_type>({ 16, 3 }) * -1)));
        }

        TEST_F(single_flow_router__raster_rook, receivers_distance)
        {
            update();

            xt::xtensor<double, 1> expected_fixed_receivers_distance{ 1.2, 1.2, 0.0, 1.1, 1.1, 1.1,
                                                                      1.1, 0.0, 1.2, 0.0, 1.2, 1.1,
                                                                      0.0, 1.1, 0.0, 1.2 };

            EXPECT_TRUE(xt::allclose(xt::col(fixed_graph.receivers_distance(), 0),
                                     expected_fixed_receivers_distance));
            EXPECT_TRUE(
                xt::all(xt::equal(xt::view(fixed_graph.receivers_distance(), 0, xt::range(1, 4)),
                                  xt::ones<double>({ 16, 3 }) * -1)));

            xt::xtensor<double, 1> expected_looped_receivers_distance{ 1.2, 1.2, 1.1, 1.1, 1.2, 1.1,
                                                                       1.1, 0.0, 1.2, 0.0, 1.2, 1.1,
                                                                       0.0, 1.1, 0.0, 1.2 };

            EXPECT_TRUE(xt::allclose(xt::col(looped_graph.receivers_distance(), 0),
                                     expected_looped_receivers_distance));
            EXPECT_TRUE(
                xt::all(xt::equal(xt::view(looped_graph.receivers_distance(), 0, xt::range(1, 4)),
                                  xt::ones<double>({ 16, 3 }) * -1)));
        }

        TEST_F(single_flow_router__raster_rook, receivers_count)
        {
            update();

            EXPECT_TRUE(
                xt::all(xt::equal(fixed_graph.receivers_count(), xt::ones<std::uint8_t>({ 16 }))));

            EXPECT_TRUE(
                xt::all(xt::equal(looped_graph.receivers_count(), xt::ones<std::uint8_t>({ 16 }))));
        }

        TEST_F(single_flow_router__raster_rook, receivers_weight)
        {
            update();

            EXPECT_TRUE(xt::all(
                xt::equal(xt::col(fixed_graph.receivers_weight(), 0), xt::ones<double>({ 16 }))));
            EXPECT_TRUE(xt::all(
                xt::equal(xt::view(fixed_graph.receivers_weight(), xt::all(), xt::range(1, 4)),
                          xt::zeros<double>({ 16, 3 }))));

            EXPECT_TRUE(xt::all(
                xt::equal(xt::col(looped_graph.receivers_weight(), 0), xt::ones<double>({ 16 }))));
            EXPECT_TRUE(xt::all(
                xt::equal(xt::view(looped_graph.receivers_weight(), xt::all(), xt::range(1, 4)),
                          xt::zeros<double>({ 16, 3 }))));
        }

        TEST_F(single_flow_router__raster_rook, donors)
        {
            update();

            xt::xarray<size_type> expected_fixed_donors = xt::ones<size_type>({ 16, 5 }) * -1;
            expected_fixed_donors(1, 0) = 0;
            xt::view(expected_fixed_donors, 2, xt::range(0, 2)) = xt::xarray<size_type>({ 1, 2 });
            xt::view(expected_fixed_donors, 7, xt::range(0, 3))
                = xt::xarray<size_type>({ 3, 7, 11 });
            expected_fixed_donors(8, 0) = 4;
            xt::view(expected_fixed_donors, 9, xt::range(0, 5))
                = xt::xarray<size_type>({ 5, 8, 9, 10, 13 });
            expected_fixed_donors(10, 0) = 6;
            expected_fixed_donors(12, 0) = 12;
            xt::view(expected_fixed_donors, 14, xt::range(0, 2))
                = xt::xarray<size_type>({ 14, 15 });

            EXPECT_TRUE(xt::all(xt::equal(fixed_graph.donors(), expected_fixed_donors)));

            xt::xarray<size_type> expected_looped_donors = xt::ones<size_type>({ 16, 5 }) * -1;
            expected_looped_donors(1, 0) = 0;
            expected_looped_donors(2, 0) = 1;
            xt::view(expected_looped_donors, 7, xt::range(0, 4))
                = xt::xarray<size_type>({ 3, 4, 7, 11 });
            xt::view(expected_looped_donors, 9, xt::range(0, 5))
                = xt::xarray<size_type>({ 5, 8, 9, 10, 13 });
            expected_looped_donors(10, 0) = 6;
            expected_looped_donors(12, 0) = 12;
            xt::view(expected_looped_donors, 14, xt::range(0, 3))
                = xt::xarray<size_type>({ 2, 14, 15 });

            EXPECT_TRUE(xt::all(xt::equal(looped_graph.donors(), expected_looped_donors)));
        }

        TEST_F(single_flow_router__raster_rook, donors_count)
        {
            update();

            xt::xtensor<std::uint8_t, 1> expected_fixed_donors_count{ 0, 1, 2, 0, 0, 0, 0, 3,
                                                                      1, 5, 1, 0, 1, 0, 2, 0 };

            EXPECT_TRUE(
                xt::all(xt::equal(fixed_graph.donors_count(), expected_fixed_donors_count)));

            xt::xtensor<std::uint8_t, 1> expected_looped_donors_count{ 0, 1, 1, 0, 0, 0, 0, 4,
                                                                       0, 5, 1, 0, 1, 0, 3, 0 };

            EXPECT_TRUE(
                xt::all(xt::equal(looped_graph.donors_count(), expected_looped_donors_count)));
        }

        TEST_F(single_flow_router__raster_rook, dfs_stack)
        {
            update();

            xt::xtensor<std::uint8_t, 1> expected_fixed_stack{ 2, 1, 0,  7, 3,  11, 9,  5,
                                                               8, 4, 10, 6, 13, 12, 14, 15 };
            EXPECT_TRUE(xt::all(xt::equal(fixed_graph.dfs_stack(), expected_fixed_stack)));

            xt::xtensor<std::uint8_t, 1> expected_looped_stack{ 7, 3,  4,  11, 9, 5, 8, 10,
                                                                6, 13, 12, 14, 2, 1, 0, 15 };
            EXPECT_TRUE(xt::all(xt::equal(looped_graph.dfs_stack(), expected_looped_stack)));
        }

        TEST_F(single_flow_router__raster_rook, dfs_iterators)
        {
            update();

            xt::xtensor<std::uint8_t, 1> expected_fixed_stack{ 2, 1, 0,  7, 3,  11, 9,  5,
                                                               8, 4, 10, 6, 13, 12, 14, 15 };

            auto bfs = fixed_graph.dfs_cbegin();
            for (std::size_t i = 0; i < fixed_graph.size(); ++i)
            {
                EXPECT_EQ(*bfs, expected_fixed_stack[i]);
                ++bfs;
            }

            auto rbfs = fixed_graph.dfs_crbegin();
            for (std::size_t i = fixed_graph.size(); i != 0; --i)
            {
                EXPECT_EQ(*rbfs, expected_fixed_stack[i - 1]);
                ++rbfs;
            }

            std::size_t sum = 0;
            for (auto it = fixed_graph.dfs_cbegin(); it != fixed_graph.dfs_cend(); ++it)
            {
                sum += *it;
            }
            EXPECT_EQ(sum, 120u);

            sum = 0;
            for (auto it = fixed_graph.dfs_crbegin(); it != fixed_graph.dfs_crend(); ++it)
            {
                sum += *it;
            }
            EXPECT_EQ(sum, 120u);
        }


        class single_flow_router__raster_bishop : public single_flow_router__raster_base
        {
        protected:
            using grid_type = fs::raster_grid_xt<fs::xtensor_selector, fs::raster_connect::bishop>;
            using flow_graph_type = fs::flow_graph<grid_type, double>;

            grid_type fixed_grid
                = grid_type(shape, { 1.1, 1.2 }, fs::node_status::fixed_value_boundary);
            grid_type looped_grid
                = grid_type(shape, { 1.1, 1.2 }, fs::node_status::looped_boundary);

            flow_graph_type fixed_graph
                = flow_graph_type(fixed_grid,
                                  std::make_unique<fs::single_flow_router<flow_graph_type>>(),
                                  std::make_unique<fs::no_sink_resolver<flow_graph_type>>());

            flow_graph_type looped_graph
                = flow_graph_type(looped_grid,
                                  std::make_unique<fs::single_flow_router<flow_graph_type>>(),
                                  std::make_unique<fs::no_sink_resolver<flow_graph_type>>());

            xt::xtensor<double, 2> elevation{ { 0.82, 0.16, 0.14, 0.20 },
                                              { 0.71, 0.97, 0.41, 0.09 },
                                              { 0.49, 0.01, 0.10, 0.38 },
                                              { 0.29, 0.82, 0.09, 0.88 } };

            void update()
            {
                fixed_graph.update_routes(elevation);
                looped_graph.update_routes(elevation);
            }
        };

        TEST_F(single_flow_router__raster_bishop, receivers)
        {
            update();

            xt::xtensor<size_type, 1> expected_fixed_receivers{ 0, 1, 7, 3,  9, 10, 9, 7,
                                                                8, 9, 7, 14, 9, 10, 9, 10 };

            EXPECT_TRUE(
                xt::all(xt::equal(xt::col(fixed_graph.receivers(), 0), expected_fixed_receivers)));
            EXPECT_TRUE(xt::all(xt::equal(xt::view(fixed_graph.receivers(), 0, xt::range(1, 4)),
                                          xt::ones<size_type>({ 16, 3 }) * -1)));

            xt::xtensor<size_type, 1> expected_looped_receivers{ 7, 14, 7, 14, 9, 10, 9, 7,
                                                                 7, 9,  7, 14, 9, 10, 9, 10 };

            EXPECT_TRUE(xt::all(
                xt::equal(xt::col(looped_graph.receivers(), 0), expected_looped_receivers)));
            EXPECT_TRUE(xt::all(xt::equal(xt::view(looped_graph.receivers(), 0, xt::range(1, 4)),
                                          xt::ones<size_type>({ 16, 3 }) * -1)));
        }

        TEST_F(single_flow_router__raster_bishop, receivers_distance)
        {
            update();

            xt::xtensor<double, 1> expected_fixed_receivers_distance{ 0.0, 0.0, dia, 0.0, dia, dia,
                                                                      dia, 0.0, 0.0, 0.0, dia, dia,
                                                                      dia, dia, dia, dia };

            EXPECT_TRUE(xt::allclose(xt::col(fixed_graph.receivers_distance(), 0),
                                     expected_fixed_receivers_distance));
            EXPECT_TRUE(
                xt::all(xt::equal(xt::view(fixed_graph.receivers_distance(), 0, xt::range(1, 4)),
                                  xt::ones<double>({ 16, 3 }) * -1)));

            xt::xtensor<double, 1> expected_looped_receivers_distance{ dia, dia, dia, dia, dia, dia,
                                                                       dia, 0.0, dia, 0.0, dia, dia,
                                                                       dia, dia, dia, dia };

            EXPECT_TRUE(xt::allclose(xt::col(looped_graph.receivers_distance(), 0),
                                     expected_looped_receivers_distance));
            EXPECT_TRUE(
                xt::all(xt::equal(xt::view(looped_graph.receivers_distance(), 0, xt::range(1, 4)),
                                  xt::ones<double>({ 16, 3 }) * -1)));
        }

        TEST_F(single_flow_router__raster_bishop, receivers_count)
        {
            update();

            EXPECT_TRUE(
                xt::all(xt::equal(fixed_graph.receivers_count(), xt::ones<std::uint8_t>({ 16 }))));

            EXPECT_TRUE(
                xt::all(xt::equal(looped_graph.receivers_count(), xt::ones<std::uint8_t>({ 16 }))));
        }

        TEST_F(single_flow_router__raster_bishop, receivers_weight)
        {
            update();

            EXPECT_TRUE(xt::all(
                xt::equal(xt::col(fixed_graph.receivers_weight(), 0), xt::ones<double>({ 16 }))));
            EXPECT_TRUE(xt::all(
                xt::equal(xt::view(fixed_graph.receivers_weight(), xt::all(), xt::range(1, 4)),
                          xt::zeros<double>({ 16, 3 }))));

            EXPECT_TRUE(xt::all(
                xt::equal(xt::col(looped_graph.receivers_weight(), 0), xt::ones<double>({ 16 }))));
            EXPECT_TRUE(xt::all(
                xt::equal(xt::view(looped_graph.receivers_weight(), xt::all(), xt::range(1, 4)),
                          xt::zeros<double>({ 16, 3 }))));
        }

        TEST_F(single_flow_router__raster_bishop, donors)
        {
            update();

            xt::xarray<size_type> expected_fixed_donors = xt::ones<size_type>({ 16, 5 }) * -1;
            expected_fixed_donors(0, 0) = 0;
            expected_fixed_donors(1, 0) = 1;
            expected_fixed_donors(3, 0) = 3;
            xt::view(expected_fixed_donors, 7, xt::range(0, 3))
                = xt::xarray<size_type>({ 2, 7, 10 });
            expected_fixed_donors(8, 0) = 8;
            xt::view(expected_fixed_donors, 9, xt::range(0, 5))
                = xt::xarray<size_type>({ 4, 6, 9, 12, 14 });
            xt::view(expected_fixed_donors, 10, xt::range(0, 3))
                = xt::xarray<size_type>({ 5, 13, 15 });
            expected_fixed_donors(14, 0) = 11;

            EXPECT_TRUE(xt::all(xt::equal(fixed_graph.donors(), expected_fixed_donors)));

            xt::xarray<size_type> expected_looped_donors = xt::ones<size_type>({ 16, 5 }) * -1;
            xt::view(expected_looped_donors, 7, xt::range(0, 5))
                = xt::xarray<size_type>({ 0, 2, 7, 8, 10 });
            xt::view(expected_looped_donors, 9, xt::range(0, 5))
                = xt::xarray<size_type>({ 4, 6, 9, 12, 14 });
            xt::view(expected_looped_donors, 10, xt::range(0, 3))
                = xt::xarray<size_type>({ 5, 13, 15 });
            xt::view(expected_looped_donors, 14, xt::range(0, 3))
                = xt::xarray<size_type>({ 1, 3, 11 });

            EXPECT_TRUE(xt::all(xt::equal(looped_graph.donors(), expected_looped_donors)));
        }

        TEST_F(single_flow_router__raster_bishop, donors_count)
        {
            update();

            xt::xtensor<std::uint8_t, 1> expected_fixed_donors_count{ 1, 1, 0, 1, 0, 0, 0, 3,
                                                                      1, 5, 3, 0, 0, 0, 1, 0 };

            EXPECT_TRUE(
                xt::all(xt::equal(fixed_graph.donors_count(), expected_fixed_donors_count)));

            xt::xtensor<std::uint8_t, 1> expected_looped_donors_count{ 0, 0, 0, 0, 0, 0, 0, 5,
                                                                       0, 5, 3, 0, 0, 0, 3, 0 };

            EXPECT_TRUE(
                xt::all(xt::equal(looped_graph.donors_count(), expected_looped_donors_count)));
        }

        TEST_F(single_flow_router__raster_bishop, dfs_stack)
        {
            update();

            xt::xtensor<std::uint8_t, 1> expected_fixed_stack{ 0,  1, 3, 7, 2, 10, 5,  13,
                                                               15, 8, 9, 4, 6, 12, 14, 11 };
            EXPECT_TRUE(xt::all(xt::equal(fixed_graph.dfs_stack(), expected_fixed_stack)));

            xt::xtensor<std::uint8_t, 1> expected_looped_stack{ 7, 0, 2, 8,  10, 5, 13, 15,
                                                                9, 4, 6, 12, 14, 1, 3,  11 };
            EXPECT_TRUE(xt::all(xt::equal(looped_graph.dfs_stack(), expected_looped_stack)));
        }

        TEST_F(single_flow_router__raster_bishop, dfs_iterators)
        {
            update();

            xt::xtensor<std::uint8_t, 1> expected_fixed_stack{ 0,  1, 3, 7, 2, 10, 5,  13,
                                                               15, 8, 9, 4, 6, 12, 14, 11 };

            auto bfs = fixed_graph.dfs_cbegin();
            for (std::size_t i = 0; i < fixed_graph.size(); ++i)
            {
                EXPECT_EQ(*bfs, expected_fixed_stack[i]);
                ++bfs;
            }

            auto rbfs = fixed_graph.dfs_crbegin();
            for (std::size_t i = fixed_graph.size(); i != 0; --i)
            {
                EXPECT_EQ(*rbfs, expected_fixed_stack[i - 1]);
                ++rbfs;
            }

            std::size_t sum = 0;
            for (auto it = fixed_graph.dfs_cbegin(); it != fixed_graph.dfs_cend(); ++it)
            {
                sum += *it;
            }
            EXPECT_EQ(sum, 120u);

            sum = 0;
            for (auto it = fixed_graph.dfs_crbegin(); it != fixed_graph.dfs_crend(); ++it)
            {
                sum += *it;
            }
            EXPECT_EQ(sum, 120u);
        }
    }
}
