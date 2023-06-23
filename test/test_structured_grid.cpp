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

            fs::node_status fixed = fs::node_status::fixed_value;
            std::array<node_s, 2> loop{ { fs::node_status::looped, fs::node_status::looped } };

            fs::profile_boundary_status fixed_status{ fixed };
            fs::profile_boundary_status looped_status{ loop };

            using grid_type = fs::profile_grid_xt<fs::xt_selector>;
            using size_type = typename grid_type::size_type;

            size_type shape{ 5 };
            grid_type fixed_grid = grid_type(shape, 1.3, fs::node_status::fixed_value);
            grid_type looped_grid = grid_type(shape, 1.4, fs::node_status::looped);
        };

        TEST_F(structured_grid, nodes_indices)
        {
            auto nodes_indices = fixed_grid.nodes_indices();

            std::size_t sum = 0;
            std::size_t size = 0;
            for (auto it = nodes_indices.begin(); it != nodes_indices.end(); ++it)
            {
                sum += *it;
                ++size;
            }
            EXPECT_EQ(sum, 10u);
            EXPECT_EQ(size, fixed_grid.size());

            sum = 0;
            size = 0;
            for (auto it = nodes_indices.rbegin(); it != nodes_indices.rend(); ++it)
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

        TEST_F(structured_grid, nodes_indices_status)
        {
            auto nodes_indices_fvalue = fixed_grid.nodes_indices(fs::node_status::fixed_value);

            {
                SCOPED_TRACE("test fixed_value begin");

                auto it = nodes_indices_fvalue.begin();

                EXPECT_EQ(*it, 0u);
                EXPECT_EQ(*(it++), 0u);
                EXPECT_EQ(*(it++), 4u);
                EXPECT_EQ(*it, 5u);
                EXPECT_EQ(*(--it), 4u);
                EXPECT_EQ(*(++it), 5u);
                EXPECT_EQ(*(++it), 6u);
            }

            {
                SCOPED_TRACE("test fixed_value end");

                auto it = nodes_indices_fvalue.end();

                EXPECT_EQ(*it, 5u);
                EXPECT_EQ(*(--it), 4u);
                EXPECT_EQ(*(--it), 0u);
            }

            auto nodes_indices_core = fixed_grid.nodes_indices(fs::node_status::core);

            {
                SCOPED_TRACE("test core begin");

                auto it = nodes_indices_core.begin();

                EXPECT_EQ(*(it++), 1u);
                EXPECT_EQ(*(it++), 2u);
                EXPECT_EQ(*it, 3u);
                EXPECT_EQ(*(--it), 2u);
            }

            {
                SCOPED_TRACE("test core end");

                auto it = nodes_indices_core.end();

                EXPECT_EQ(*it, 5u);
                EXPECT_EQ(*(--it), 3u);
            }
        }

    }
}
