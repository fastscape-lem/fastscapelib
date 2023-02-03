#include "fastscapelib/grid/structured_grid.hpp"
#include "fastscapelib/grid/profile_grid.hpp"

#include "gtest/gtest.h"


namespace fs = fastscapelib;


namespace fastscapelib
{
    namespace testing
    {

        class neighbor : public ::testing::Test
        {
        protected:
            fs::neighbor n{ 3, 1.35, fs::node_status::core };
        };

        TEST_F(neighbor, ctor)
        {
            EXPECT_EQ(n.idx, 3u);
            EXPECT_EQ(n.distance, 1.35);
            EXPECT_EQ(n.status, fs::node_status::core);
        }

        TEST_F(neighbor, equal)
        {
            fs::neighbor other_n{ 3, 1.35, fs::node_status::core };
            EXPECT_EQ(n, other_n);
        }


        class structured_grid : public ::testing::Test
        {
        protected:
            using node_s = fs::node_status;

            fs::node_status fixed = fs::node_status::fixed_value_boundary;
            std::array<node_s, 2> loop{ { fs::node_status::looped_boundary,
                                          fs::node_status::looped_boundary } };

            fs::profile_boundary_status fixed_status{ fixed };
            fs::profile_boundary_status looped_status{ loop };

            using grid_type = fs::profile_grid_xt<fs::xt_selector>;
            using size_type = typename grid_type::size_type;

            size_type shape{ 5 };
            grid_type fixed_grid = grid_type(shape, 1.3, fs::node_status::fixed_value_boundary);
            grid_type looped_grid = grid_type(shape, 1.4, fs::node_status::looped_boundary);
        };

        TEST_F(structured_grid, index_iterator)
        {
            index_iterator<grid_type> it(fixed_grid);

            EXPECT_EQ(*(it++), 0u);
            EXPECT_EQ(*(it++), 1u);
        }

        TEST_F(structured_grid, nodes_indices_begin)
        {
            auto boundaries_filter
                = fixed_grid.nodes_indices_begin(fs::node_status::fixed_value_boundary);
            EXPECT_EQ(*boundaries_filter, 0u);
            EXPECT_EQ(*(boundaries_filter++), 0u);
            EXPECT_EQ(*(boundaries_filter++), 4u);
            EXPECT_EQ(*boundaries_filter, 5u);
            EXPECT_EQ(*(--boundaries_filter), 4u);
            EXPECT_EQ(*(++boundaries_filter), 5u);
            EXPECT_EQ(*(++boundaries_filter), 6u);

            auto core_filter = fixed_grid.nodes_indices_begin(fs::node_status::core);
            EXPECT_EQ(*(core_filter++), 1u);
            EXPECT_EQ(*(core_filter++), 2u);
            EXPECT_EQ(*core_filter, 3u);
            EXPECT_EQ(*(--core_filter), 2u);
        }

        TEST_F(structured_grid, nodes_indices_end)
        {
            auto boundaries_filter
                = fixed_grid.nodes_indices_end(fs::node_status::fixed_value_boundary);
            EXPECT_EQ(*boundaries_filter, 5u);
            EXPECT_EQ(*(--boundaries_filter), 4u);
            EXPECT_EQ(*(--boundaries_filter), 0u);

            auto core_filter_end = fixed_grid.nodes_indices_end(fs::node_status::core);
            EXPECT_EQ(*core_filter_end, 5u);
            EXPECT_EQ(*(--core_filter_end), 3u);
        }

        TEST_F(structured_grid, nodes_indices)
        {
            auto node_indices = fixed_grid.nodes_indices();

            std::size_t sum = 0;
            std::size_t size = 0;
            for (auto it = node_indices.begin(); it != node_indices.end(); ++it)
            {
                sum += *it;
                ++size;
            }
            EXPECT_EQ(sum, 10u);
            EXPECT_EQ(size, fixed_grid.size());

            sum = 0;
            size = 0;
            for (auto it = node_indices.rbegin(); it != node_indices.rend(); ++it)
            {
                sum += *it;
                ++size;
            }
            EXPECT_EQ(sum, 10u);
            EXPECT_EQ(size, fixed_grid.size());

            sum = 0;
            size = 0;
            for (auto idx : fixed_grid.nodes_indices())
            {
                sum += idx;
                ++size;
            }
            EXPECT_EQ(sum, 10u);
            EXPECT_EQ(size, fixed_grid.size());
        }
    }
}
