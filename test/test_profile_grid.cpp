#include "fastscapelib/profile_grid.hpp"

#include "gtest/gtest.h"


namespace fs = fastscapelib;


namespace fastscapelib
{
    namespace testing
    {

        class neighbor: public ::testing::Test
        {
            protected:

                fs::neighbor n {3, 1.35, fs::node_status::core};
        };

        TEST_F(neighbor, ctor)
        {
            EXPECT_EQ(n.idx, 3u);
            EXPECT_EQ(n.distance, 1.35);
            EXPECT_EQ(n.status, fs::node_status::core);
        }

        TEST_F(neighbor, equal)
        {
            fs::neighbor other_n {3, 1.35, fs::node_status::core};
            EXPECT_EQ(n, other_n);
        }


        class profile_boundary_status: public ::testing::Test
        {
            protected:

                using node_s = fs::node_status;

                node_s fixed = node_s::fixed_value_boundary;
                std::array<node_s, 2> loop {{node_s::looped_boundary, node_s::looped_boundary}};
                std::array<node_s, 2> hill_formed_loop {{node_s::looped_boundary, node_s::fixed_value_boundary}};

                fs::profile_boundary_status fixed_value_status {fixed};
                fs::profile_boundary_status looped_status {loop};
                
        };

        TEST_F(profile_boundary_status, ctor)
        {
            EXPECT_EQ(fixed_value_status.left, fixed);
            EXPECT_EQ(fixed_value_status.right, fixed);

            EXPECT_EQ(looped_status.left, node_s::looped_boundary);
            EXPECT_EQ(looped_status.right, node_s::looped_boundary);

            EXPECT_THROW(fs::profile_boundary_status {hill_formed_loop}, std::invalid_argument);
        }

        TEST_F(profile_boundary_status, is_horizontal_looped)
        {
            EXPECT_FALSE(fixed_value_status.is_horizontal_looped());
            EXPECT_TRUE(looped_status.is_horizontal_looped());
        }


        class profile_grid: public ::testing::Test
        {
            protected:

                using node_s = fs::node_status;

                node_s fixed = node_s::fixed_value_boundary;
                std::array<node_s, 2> loop {{node_s::looped_boundary, node_s::looped_boundary}};
                std::array<node_s, 2> hill_formed_loop {{node_s::looped_boundary, node_s::fixed_value_boundary}};

                fs::profile_boundary_status fixed_value_status {fixed};
                fs::profile_boundary_status looped_status {loop};

                using grid_type = fs::profile_grid_xt<fs::xtensor_selector>;
                using size_type = typename grid_type::size_type;

                size_type shape {5};
                grid_type fixed_grid = grid_type(shape, 1.3, fs::node_status::fixed_value_boundary);
                grid_type looped_grid = grid_type(shape, 1.4, fs::node_status::looped_boundary);

                fs::node_status status_fixed(std::size_t idx)
                { 
                    return ((idx == 0) || (idx == 4)) ? fs::node_status::fixed_value_boundary : fs::node_status::core; 
                };

                fs::node_status status_looped(std::size_t idx)
                { 
                    return ((idx == 0) || (idx == 4)) ? fs::node_status::looped_boundary : fs::node_status::core; 
                };
        };

        TEST_F(profile_grid, boundary_status)
        {
            fs::profile_boundary_status looped_status {node_s::looped_boundary, node_s::looped_boundary};
            ASSERT_THROW(fs::profile_boundary_status status {hill_formed_loop}, std::invalid_argument);

            EXPECT_FALSE(fixed_value_status.is_horizontal_looped());
            EXPECT_TRUE(looped_status.is_horizontal_looped());
        }

        TEST_F(profile_grid, ctor)
        {
            std::vector<fs::node> nodes_vector1 {fs::node({1, fs::node_status::fixed_value_boundary})};
            grid_type g1(shape, 1.3, fs::node_status::looped_boundary, nodes_vector1);

            auto expected_status = grid_type::node_status_type {fs::node_status::looped_boundary, 
                                                                fs::node_status::fixed_value_boundary, 
                                                                fs::node_status::core, 
                                                                fs::node_status::core, 
                                                                fs::node_status::looped_boundary};
            ASSERT_EQ(g1.status_at_nodes(), expected_status);

            std::vector<fs::node> nodes_vector2 {fs::node({15, fs::node_status::core})};
            ASSERT_THROW(grid_type(shape, 1.3, fs::node_status::fixed_value_boundary, nodes_vector2), std::out_of_range);

            std::vector<fs::node> nodes_vector3 {fs::node({0, fs::node_status::core})};
            ASSERT_THROW(grid_type(shape, 1.3, fs::node_status::looped_boundary, nodes_vector3), std::invalid_argument);
        }

        TEST_F(profile_grid, neighbors__fixed_value_boundary)
        {
            EXPECT_EQ(fixed_grid.neighbors_indices_cache().cache_used(), 0u);
            EXPECT_EQ(fixed_grid.neighbors_indices_cache().cache_size(), 5u);

            EXPECT_EQ(fixed_grid.neighbors(0), 
                     (xt::xtensor<fs::neighbor, 1> { {1, 1.3, fs::node_status::core} } ));

            for(std::size_t i=1; i<4; ++i)
            {
                EXPECT_EQ(fixed_grid.neighbors_indices_cache().cache_used(), i);
                EXPECT_EQ(fixed_grid.neighbors(i),
                          (xt::xtensor<fs::neighbor, 1> { {i-1, 1.3, status_fixed(i-1)},
                                                          {i+1, 1.3, status_fixed(i+1)} } ));
            }
            EXPECT_EQ(fixed_grid.neighbors_indices_cache().cache_used(), 4u);

            EXPECT_EQ(fixed_grid.neighbors(4),
                      (xt::xtensor<fs::neighbor, 1> { {3, 1.3, fs::node_status::core} } ));
            EXPECT_EQ(fixed_grid.neighbors_indices_cache().cache_used(), 5u);
        }

        TEST_F(profile_grid, neighbors__looped_boundary)
        {
            EXPECT_EQ(looped_grid.neighbors_indices_cache().cache_used(), 0u);
            EXPECT_EQ(looped_grid.neighbors_indices_cache().cache_size(), 5u);

            EXPECT_EQ(looped_grid.neighbors(0), (xt::xtensor<fs::neighbor, 1> {{4, 1.4, fs::node_status::looped_boundary},
                                                                               {1, 1.4, fs::node_status::core}}));

            for(std::size_t i=1; i<4; ++i)
            {
                EXPECT_EQ(looped_grid.neighbors_indices_cache().cache_used(), i);
                EXPECT_EQ(looped_grid.neighbors(i), (xt::xtensor<fs::neighbor, 1> {{i-1, 1.4, status_looped(i-1)},
                                                                                   {i+1, 1.4, status_looped(i+1)}}));
            }
            EXPECT_EQ(looped_grid.neighbors_indices_cache().cache_used(), 4u);
            EXPECT_EQ(looped_grid.neighbors(4), (xt::xtensor<fs::neighbor, 1> {{3, 1.4, fs::node_status::core},
                                                                               {0, 1.4, fs::node_status::looped_boundary}}));
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

        TEST_F(profile_grid, length)
        {
            EXPECT_EQ(fixed_grid.length(), 5.2);
            EXPECT_EQ(looped_grid.length(), 5.6);
        }

        TEST_F(profile_grid, from_length)
        {
            auto grid_from_length = grid_type::from_length(151, 1500., fs::node_status::fixed_value_boundary);
            EXPECT_EQ(grid_from_length.length(), 1500.);
            EXPECT_EQ(grid_from_length.size(), 151u);
            EXPECT_EQ(grid_from_length.spacing(), 10.);
        }
    }
}