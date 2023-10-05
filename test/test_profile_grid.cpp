#include <algorithm>
#include <map>
#include <vector>

#include "fastscapelib/grid/profile_grid.hpp"

#include "gtest/gtest.h"


namespace fs = fastscapelib;


namespace fastscapelib
{
    namespace testing
    {

        class profile_boundary_status : public ::testing::Test
        {
        protected:
            using node_s = fs::node_status;

            node_s fixed = node_s::fixed_value;
            std::array<node_s, 2> loop{ { node_s::looped, node_s::looped } };
            std::array<node_s, 2> hill_formed_loop{ { node_s::looped, node_s::fixed_value } };

            fs::profile_boundary_status fixed_value_status{ fixed };
            fs::profile_boundary_status looped_status{ loop };
        };

        TEST_F(profile_boundary_status, ctor)
        {
            EXPECT_EQ(fixed_value_status.left, fixed);
            EXPECT_EQ(fixed_value_status.right, fixed);

            EXPECT_EQ(looped_status.left, node_s::looped);
            EXPECT_EQ(looped_status.right, node_s::looped);

            EXPECT_THROW(fs::profile_boundary_status{ hill_formed_loop }, std::invalid_argument);
        }

        TEST_F(profile_boundary_status, is_horizontal_looped)
        {
            EXPECT_FALSE(fixed_value_status.is_horizontal_looped());
            EXPECT_TRUE(looped_status.is_horizontal_looped());
        }


        class profile_grid : public ::testing::Test
        {
        protected:
            using node_s = fs::node_status;

            node_s fixed = node_s::fixed_value;
            std::array<node_s, 2> loop{ { node_s::looped, node_s::looped } };
            std::array<node_s, 2> hill_formed_loop{ { node_s::looped, node_s::fixed_value } };

            fs::profile_boundary_status fixed_value_status{ fixed };
            fs::profile_boundary_status looped_status{ loop };

            using grid_type = fs::profile_grid_xt<fs::xt_selector>;
            using size_type = typename grid_type::size_type;
            using shape_type = typename grid_type::shape_type;

            size_type size = 5;
            grid_type fixed_grid = grid_type(size, 1.3, fs::node_status::fixed_value);
            grid_type looped_grid = grid_type(size, 1.4, fs::node_status::looped);

            fs::node_status status_fixed(std::size_t idx)
            {
                return ((idx == 0) || (idx == 4)) ? fs::node_status::fixed_value
                                                  : fs::node_status::core;
            };

            fs::node_status status_looped(std::size_t idx)
            {
                return ((idx == 0) || (idx == 4)) ? fs::node_status::looped : fs::node_status::core;
            };
        };

        TEST_F(profile_grid, boundary_status)
        {
            fs::profile_boundary_status looped_status{ node_s::looped, node_s::looped };
            ASSERT_THROW(fs::profile_boundary_status status{ hill_formed_loop },
                         std::invalid_argument);

            EXPECT_FALSE(fixed_value_status.is_horizontal_looped());
            EXPECT_TRUE(looped_status.is_horizontal_looped());
        }

        TEST_F(profile_grid, static_expr)
        {
            EXPECT_EQ(fs::profile_grid::is_structured(), true);
            EXPECT_EQ(fs::profile_grid::is_uniform(), true);
            EXPECT_EQ(fs::profile_grid::n_neighbors_max(), 2u);
            EXPECT_EQ(fs::profile_grid::xt_ndims(), 1);
        }

        TEST_F(profile_grid, ctor)
        {
            std::map<size_type, fs::node_status> nodes_status1{ { 15, fs::node_status::core } };
            ASSERT_THROW(grid_type(size, 1.3, fs::node_status::fixed_value, nodes_status1),
                         std::out_of_range);

            std::map<size_type, fs::node_status> nodes_status2{ { 0, fs::node_status::core } };
            ASSERT_THROW(grid_type(size, 1.3, fs::node_status::looped, nodes_status2),
                         std::invalid_argument);
        }

        TEST_F(profile_grid, neighbors__fixed_value)
        {
            EXPECT_EQ(fixed_grid.neighbors_indices_cache().cache_used(), 0u);
            EXPECT_EQ(fixed_grid.neighbors_indices_cache().cache_size(), 5u);

            EXPECT_EQ(fixed_grid.neighbors(0),
                      (std::vector<fs::neighbor>{ { 1, 1.3, fs::node_status::core } }));

            for (std::size_t i = 1; i < 4; ++i)
            {
                EXPECT_EQ(fixed_grid.neighbors_indices_cache().cache_used(), i);
                EXPECT_EQ(fixed_grid.neighbors(i),
                          (std::vector<fs::neighbor>{ { i - 1, 1.3, status_fixed(i - 1) },
                                                      { i + 1, 1.3, status_fixed(i + 1) } }));
            }
            EXPECT_EQ(fixed_grid.neighbors_indices_cache().cache_used(), 4u);

            EXPECT_EQ(fixed_grid.neighbors(4),
                      (std::vector<fs::neighbor>{ { 3, 1.3, fs::node_status::core } }));
            EXPECT_EQ(fixed_grid.neighbors_indices_cache().cache_used(), 5u);
        }

        TEST_F(profile_grid, neighbors__looped)
        {
            EXPECT_EQ(looped_grid.neighbors_indices_cache().cache_used(), 0u);
            EXPECT_EQ(looped_grid.neighbors_indices_cache().cache_size(), 5u);

            EXPECT_EQ(looped_grid.neighbors(0),
                      (std::vector<fs::neighbor>{ { 4, 1.4, fs::node_status::looped },
                                                  { 1, 1.4, fs::node_status::core } }));

            for (std::size_t i = 1; i < 4; ++i)
            {
                EXPECT_EQ(looped_grid.neighbors_indices_cache().cache_used(), i);
                EXPECT_EQ(looped_grid.neighbors(i),
                          (std::vector<fs::neighbor>{ { i - 1, 1.4, status_looped(i - 1) },
                                                      { i + 1, 1.4, status_looped(i + 1) } }));
            }
            EXPECT_EQ(looped_grid.neighbors_indices_cache().cache_used(), 4u);
            EXPECT_EQ(looped_grid.neighbors(4),
                      (std::vector<fs::neighbor>{ { 3, 1.4, fs::node_status::core },
                                                  { 0, 1.4, fs::node_status::looped } }));
            EXPECT_EQ(looped_grid.neighbors_indices_cache().cache_used(), 5u);
        }

        TEST_F(profile_grid, spacing)
        {
            EXPECT_EQ(fixed_grid.spacing(), 1.3);
            EXPECT_EQ(looped_grid.spacing(), 1.4);
        }

        TEST_F(profile_grid, size)
        {
            EXPECT_EQ(fixed_grid.size(), 5u);
            EXPECT_EQ(looped_grid.size(), 5u);
        }

        TEST_F(profile_grid, shape)
        {
            EXPECT_EQ(fixed_grid.shape(), shape_type({ 5 }));
            EXPECT_EQ(looped_grid.shape(), shape_type({ 5 }));
        }

        TEST_F(profile_grid, length)
        {
            EXPECT_EQ(fixed_grid.length(), 5.2);
            EXPECT_EQ(looped_grid.length(), 5.6);
        }

        TEST_F(profile_grid, nodes_indices)
        {
            {
                SCOPED_TRACE("no filter");
                auto indices = fixed_grid.nodes_indices();
                std::vector<size_type> actual;
                std::copy(indices.begin(), indices.end(), std::back_inserter(actual));
                std::vector<size_type> expected{ 0, 1, 2, 3, 4 };
                EXPECT_EQ(actual, expected);
            }
            {
                SCOPED_TRACE("filter (not empty result)");
                auto indices = fixed_grid.nodes_indices(node_status::fixed_value);
                std::vector<size_type> actual;
                std::copy(indices.begin(), indices.end(), std::back_inserter(actual));
                std::vector<size_type> expected{ 0, 4 };
                EXPECT_EQ(actual, expected);
            }
            {
                SCOPED_TRACE("filter (empty result)");
                auto indices = fixed_grid.nodes_indices(node_status::fixed_gradient);
                std::vector<size_type> actual;
                std::copy(indices.begin(), indices.end(), std::back_inserter(actual));
                EXPECT_EQ(actual.size(), 0);
            }
        }

        TEST_F(profile_grid, nodes_status)
        {
            std::map<size_type, fs::node_status> nodes_status{ { 1,
                                                                 fs::node_status::fixed_value } };
            grid_type g1(size, 1.3, fs::node_status::looped, nodes_status);

            auto expected_status = grid_type::nodes_status_type{ fs::node_status::looped,
                                                                 fs::node_status::fixed_value,
                                                                 fs::node_status::core,
                                                                 fs::node_status::core,
                                                                 fs::node_status::looped };
            ASSERT_EQ(g1.nodes_status(), expected_status);
            ASSERT_EQ(g1.nodes_status(0), fs::node_status::looped);
            ASSERT_EQ(g1.nodes_status(1), fs::node_status::fixed_value);
        }

        TEST_F(profile_grid, nodes_areas)
        {
            for (auto n : fixed_grid.nodes_indices())
            {
                EXPECT_EQ(fixed_grid.nodes_areas(n), 1.3);
                EXPECT_EQ(looped_grid.nodes_areas(n), 1.4);
            }

            EXPECT_EQ(fixed_grid.nodes_areas(), xt::ones<double>(fixed_grid.shape()) * 1.3);
            EXPECT_EQ(looped_grid.nodes_areas(), xt::ones<double>(looped_grid.shape()) * 1.4);
        }

        TEST_F(profile_grid, from_length)
        {
            auto grid_from_length
                = grid_type::from_length(151, 1500., fs::node_status::fixed_value);
            EXPECT_EQ(grid_from_length.length(), 1500.);
            EXPECT_EQ(grid_from_length.size(), 151u);
            EXPECT_EQ(grid_from_length.spacing(), 10.);
        }
    }
}
