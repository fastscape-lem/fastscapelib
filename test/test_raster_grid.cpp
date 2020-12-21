#include "fastscapelib/raster_grid.hpp"

#include "gtest/gtest.h"

#include <array>


namespace fs = fastscapelib;


namespace fastscapelib
{
    namespace testing
    {

        class raster_neighbor: public ::testing::Test
        {
            protected:

                fs::raster_neighbor n {3, 1, 6, 1.35, fs::node_status::core};
        };

        TEST_F(raster_neighbor, ctor)
        {
            EXPECT_EQ(n.flatten_idx, 3);
            EXPECT_EQ(n.row, 1);
            EXPECT_EQ(n.col, 6);
            EXPECT_EQ(n.distance, 1.35);
            EXPECT_EQ(n.status, fs::node_status::core);
        }

        TEST_F(raster_neighbor, equal)
        {
            fs::raster_neighbor other_n {3, 1, 6, 1.35, fs::node_status::core};
            EXPECT_EQ(n, other_n);
        }

       class raster_boundary_status: public ::testing::Test
        {
            protected:

                using node_s = fs::node_status;

                node_s fixed = node_s::fixed_value_boundary;
                std::array<node_s, 4> hloop {{node_s::looped_boundary, node_s::looped_boundary, node_s::fixed_value_boundary, node_s::fixed_value_boundary}};
                std::array<node_s, 4> vloop {{node_s::fixed_value_boundary, node_s::fixed_value_boundary, node_s::looped_boundary, node_s::looped_boundary}};
                std::array<node_s, 4> hvloop {{node_s::looped_boundary, node_s::looped_boundary, node_s::looped_boundary, node_s::looped_boundary}};
                std::array<node_s, 4> hill_formed_loop {{node_s::looped_boundary, node_s::fixed_value_boundary, node_s::looped_boundary, node_s::looped_boundary}};

                fs::raster_boundary_status fixed_value_status {fixed};
                fs::raster_boundary_status hlooped_status {hloop};
                fs::raster_boundary_status vlooped_status {vloop};
                fs::raster_boundary_status hvlooped_status {hvloop};
        };

        TEST_F(raster_boundary_status, ctor)
        {
            EXPECT_EQ(fixed_value_status.left, fixed);
            EXPECT_EQ(fixed_value_status.right, fixed);
            EXPECT_EQ(fixed_value_status.bottom, fixed);
            EXPECT_EQ(fixed_value_status.top, fixed);

            EXPECT_EQ(hlooped_status.left, node_s::looped_boundary);
            EXPECT_EQ(hlooped_status.right, node_s::looped_boundary);
            EXPECT_EQ(hlooped_status.bottom, node_s::fixed_value_boundary);
            EXPECT_EQ(hlooped_status.top, node_s::fixed_value_boundary);

            EXPECT_EQ(vlooped_status.left, node_s::fixed_value_boundary);
            EXPECT_EQ(vlooped_status.right, node_s::fixed_value_boundary);
            EXPECT_EQ(vlooped_status.bottom, node_s::looped_boundary);
            EXPECT_EQ(vlooped_status.top, node_s::looped_boundary);

            EXPECT_EQ(hvlooped_status.left, node_s::looped_boundary);
            EXPECT_EQ(hvlooped_status.right, node_s::looped_boundary);
            EXPECT_EQ(hvlooped_status.bottom, node_s::looped_boundary);
            EXPECT_EQ(hvlooped_status.top, node_s::looped_boundary);

            EXPECT_THROW(fs::raster_boundary_status {hill_formed_loop}, std::invalid_argument);
        }

        TEST_F(raster_boundary_status, is_horizontal_looped)
        {
            EXPECT_FALSE(fixed_value_status.is_horizontal_looped());
            EXPECT_TRUE(hlooped_status.is_horizontal_looped());
            EXPECT_FALSE(vlooped_status.is_horizontal_looped());
            EXPECT_TRUE(hvlooped_status.is_horizontal_looped());
        }

        TEST_F(raster_boundary_status, is_vertical_looped)
        {
            EXPECT_FALSE(fixed_value_status.is_vertical_looped());
            EXPECT_FALSE(hlooped_status.is_vertical_looped());
            EXPECT_TRUE(vlooped_status.is_vertical_looped());
            EXPECT_TRUE(hvlooped_status.is_vertical_looped());
        }

        class raster_grid: public ::testing::Test
        {
            protected:
                using node_s = fs::node_status;

                node_s fixed = node_s::fixed_value_boundary;
                std::array<node_s, 4> hloop {{ node_s::looped_boundary, node_s::looped_boundary, node_s::fixed_value_boundary, node_s::fixed_value_boundary }};
                std::array<node_s, 4> vloop {{ node_s::fixed_value_boundary, node_s::fixed_value_boundary, node_s::looped_boundary, node_s::looped_boundary }};
                std::array<node_s, 4> hvloop {{ node_s::looped_boundary, node_s::looped_boundary, node_s::looped_boundary, node_s::looped_boundary }};
                std::array<node_s, 4> hill_formed_loop {{ node_s::looped_boundary, node_s::fixed_value_boundary, node_s::looped_boundary, node_s::looped_boundary }};

                fs::raster_boundary_status fixed_value_status {fixed};
                fs::raster_boundary_status hlooped_status {hloop};
                fs::raster_boundary_status vlooped_status {vloop};
                fs::raster_boundary_status hvlooped_status {hvloop};

                using grid_type = fs::raster_grid_xt<fs::xtensor_selector>;
                using size_type = typename grid_type::size_type;

                std::array<size_type, 2> shape {{5, 10}};
                grid_type fixed_grid = grid_type(shape, {1.3, 1.2}, fs::node_status::fixed_value_boundary);
                grid_type looped_grid = grid_type(shape, {1.4, 1.8}, fs::node_status::looped_boundary);

                fs::node_status status_fixed(std::size_t row, std::size_t col)
                { 
                    return ((row == 0) || (row == 4) || (col == 0) || (col == 9)) ? fs::node_status::fixed_value_boundary : fs::node_status::core; 
                };

                fs::node_status status_looped(std::size_t row, std::size_t col)
                { 
                    return ((row == 0) || (row == 4) || (col == 0) || (col == 9)) ? fs::node_status::looped_boundary : fs::node_status::core; 
                };
        };

        TEST_F(raster_grid, ctor)
        {
            std::array<size_type, 2> shape {{ 3, 3 }};
            std::vector<fs::raster_node> nodes_vector1 {fs::raster_node({1, 1, fs::node_status::fixed_gradient_boundary})};
            grid_type g1 = grid_type(shape, {1.4, 1.8}, hloop, nodes_vector1);

            auto expected_status = grid_type::node_status_type {{fs::node_status::fixed_value_boundary, fs::node_status::fixed_value_boundary, fs::node_status::fixed_value_boundary},
                                                                {fs::node_status::looped_boundary, fs::node_status::fixed_gradient_boundary, fs::node_status::looped_boundary},
                                                                {fs::node_status::fixed_value_boundary, fs::node_status::fixed_value_boundary, fs::node_status::fixed_value_boundary}};
            ASSERT_EQ(g1.status_at_nodes(), expected_status);

            std::vector<fs::raster_node> nodes_vector2 {fs::raster_node({15, 15, fs::node_status::core})};
            ASSERT_THROW(grid_type(shape, {1.4, 1.8}, fixed, nodes_vector2), std::out_of_range);

            std::vector<fs::raster_node> nodes_vector3 {fs::raster_node({0, 0, fs::node_status::core})};
            ASSERT_THROW(grid_type(shape, {1.4, 1.8}, hvloop, nodes_vector3), std::invalid_argument);
        }

        TEST_F(raster_grid, neighbor_indices__fixed_value_boundary)
        {
            { // First row
                // Top-left corner
                EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::queen>(0, 0)), (grid_type::neighbor_indices_type { /*(r,c)*/ { 0, 1 }, 
                                                                                                                              { 1, 0 }, { 1, 1 } }));
                EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::rook>(0, 0)), (grid_type::neighbor_indices_type { /*(r,c)*/ { 0, 1 }, 
                                                                                                                              { 1, 0 }          }));                                                                                                                              
                // Inner cols
                for(std::size_t c=1; c<9; ++c)
                {
                    EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::queen>(0, c)), (grid_type::neighbor_indices_type { { 0, c-1 }, /*(r,c)*/ { 0, c+1 }, 
                                                                                                                                  { 1, c-1 }, { 1, c }, { 1, c+1 } }));
                    EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::rook>(0, c)), (grid_type::neighbor_indices_type { { 0, c-1 }, /*(r,c)*/ { 0, c+1 }, 
                                                                                                                                             { 1, c },            }));
                }

                // Top-right corner
                EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::queen>(0, 9)), (grid_type::neighbor_indices_type { { 0, 8 }, /*(r,c)*/
                                                                                                                              { 1, 8 }, { 1, 9 } }));
                EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::rook>(0, 9)), (grid_type::neighbor_indices_type { { 0, 8 }, /*(r,c)*/
                                                                                                                                       { 1, 9 } }));
            }
            { // Inners rows
                for(std::size_t r=1; r<4; ++r)
                {
                    // First col
                    EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::queen>(r, 0)), (grid_type::neighbor_indices_type { { r-1, 0 }, { r-1, 1 },
                                                                                                                                   /*(r,c)*/  { r,   1 },
                                                                                                                                  { r+1, 0 }, { r+1, 1 } }));
                    EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::rook>(r, 0)), (grid_type::neighbor_indices_type { { r-1, 0 },
                                                                                                                                  /*(r,c)*/  { r,   1 },
                                                                                                                                 { r+1, 0 }             }));                                                                                                                                  
                    // Inners cols
                    for(std::size_t c=1; c<9; ++c)
                    {
                        EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::queen>(r, c)), (grid_type::neighbor_indices_type { { r-1, c-1 }, { r-1, c }, { r-1, c+1 },
                                                                                                                                      { r,   c-1 },  /*(r,c)*/  { r,   c+1 },
                                                                                                                                      { r+1, c-1 }, { r+1, c }, { r+1, c+1 } }));
                        EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::rook>(r, c)), (grid_type::neighbor_indices_type {               { r-1, c },
                                                                                                                                     { r,   c-1 },  /*(r,c)*/  { r,   c+1 },
                                                                                                                                                   { r+1, c },              }));
                    }

                    // last col
                    EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::queen>(r, 9)), (grid_type::neighbor_indices_type { { r-1, 8 }, { r-1, 9 },
                                                                                                                                  { r,   8 },  /*(r,c)*/
                                                                                                                                  { r+1, 8 }, { r+1, 9 } }));
                    EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::rook>(r, 9)), (grid_type::neighbor_indices_type {             { r-1, 9 }, 
                                                                                                                                 { r,   8 },  /*(r,c)*/
                                                                                                                                             { r+1, 9 } }));
                }
            }

            { // Last row
                // Bottom-left corner
                EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::queen>(4, 0)), (grid_type::neighbor_indices_type { { 3, 0 }, { 3, 1 },
                                                                                                                              /*(r,c)*/ { 4, 1 } }));
                EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::rook>(4, 0)), (grid_type::neighbor_indices_type { { 3, 0 },
                                                                                                                             /*(r,c)*/ { 4, 1 } }));
                // Inner cols
                for(std::size_t c=1; c<9; ++c)
                {
                    EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::queen>(4, c)), (grid_type::neighbor_indices_type { { 3, c-1 }, { 3, c }, { 3, c+1 },
                                                                                                                                  { 4, c-1 }, /*(r,c)*/ { 4, c+1 } }));
                    EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::rook>(4, c)), (grid_type::neighbor_indices_type {             { 3, c },
                                                                                                                                 { 4, c-1 }, /*(r,c)*/ { 4, c+1 } }));                                                                                                                                  
                }

                // Bottom-right corner
                EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::queen>(4, 9)), (grid_type::neighbor_indices_type { { 3, 8 }, { 3, 9 },
                                                                                                                              { 4, 8 }  /*(r,c)*/ }));
                EXPECT_EQ((fixed_grid.neighbor_indices<fs::raster_connect::rook>(4, 9)), (grid_type::neighbor_indices_type {          { 3, 9 },
                                                                                                                             { 4, 8 } /*(r,c)*/ }));
            }
        }

        TEST_F(raster_grid, neighbor_indices__looped_boundary)
        {
            { // First row
                // Top-left corner
                EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::queen>(0, 0)), (grid_type::neighbor_indices_type { { 4, 9 }, { 4, 0 }, { 4, 1 },
                                                                                                                               { 0, 9 }, /*(r,c)*/ { 0, 1 },
                                                                                                                               { 1, 9 }, { 1, 0 }, { 1, 1 } }));
                EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::rook>(0, 0)), (grid_type::neighbor_indices_type {           { 4, 0 },
                                                                                                                              { 0, 9 }, /*(r,c)*/ { 0, 1 },
                                                                                                                                        { 1, 0 }            }));
                // Inner cols
                for(std::size_t c=1; c<9; ++c)
                {
                    EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::queen>(0, c)), (grid_type::neighbor_indices_type { { 4, c-1 }, { 4, c }, { 4, c+1 },
                                                                                                                                   { 0, c-1 }, /*(r,c)*/ { 0, c+1 },
                                                                                                                                   { 1, c-1 }, { 1, c }, { 1, c+1 } }));
                    EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::rook>(0, c)), (grid_type::neighbor_indices_type {             { 4, c },
                                                                                                                                  { 0, c-1 }, /*(r,c)*/ { 0, c+1 },
                                                                                                                                              { 1, c }              }));
                }
                // Top-right corner
                EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::queen>(0, 9)), (grid_type::neighbor_indices_type { { 4, 8 }, { 4, 9 }, { 4, 0 },
                                                                                                                               { 0, 8 }, /*(r,c)*/ { 0, 0 },
                                                                                                                               { 1, 8 }, { 1, 9 }, { 1, 0 } }));
                EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::rook>(0, 9)), (grid_type::neighbor_indices_type {            { 4, 9 },
                                                                                                                               { 0, 8 }, /*(r,c)*/ { 0, 0 },
                                                                                                                                         { 1, 9 }            }));
            }
            { // Inners rows
                for(std::size_t r=1; r<4; ++r)
                {
                    // First col
                    EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::queen>(r, 0)), (grid_type::neighbor_indices_type { { r-1, 9 }, { r-1, 0 }, { r-1, 1 },
                                                                                                                                   { r,   9 },  /*(r,c)*/  { r,   1 },
                                                                                                                                   { r+1, 9 }, { r+1, 0 }, { r+1, 1 } }));
                    EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::rook>(r, 0)), (grid_type::neighbor_indices_type {             { r-1, 0 },
                                                                                                                                  { r,   9 },  /*(r,c)*/  { r,   1 },
                                                                                                                                              { r+1, 0 }             }));
                    // Inners cols
                    for(std::size_t c=1; c<9; ++c)
                    {
                        EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::queen>(r, c)), (grid_type::neighbor_indices_type { { r-1, c-1 }, { r-1, c }, { r-1, c+1 },
                                                                                                                                       { r,   c-1 },  /*(r,c)*/  { r,   c+1 },
                                                                                                                                       { r+1, c-1 }, { r+1, c }, { r+1, c+1 } }));
                        EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::rook>(r, c)), (grid_type::neighbor_indices_type {               { r-1, c },
                                                                                                                                      { r,   c-1 },  /*(r,c)*/ { r,   c+1 },
                                                                                                                                                    { r+1, c }               }));
                    }
                    // last col
                    EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::queen>(r, 9)), (grid_type::neighbor_indices_type { { r-1, 8 }, { r-1, 9 }, { r-1, 0 },
                                                                                                                                   { r,   8 },  /*(r,c)*/ { r,   0 },
                                                                                                                                   { r+1, 8 }, { r+1, 9 }, { r+1, 0 } }));
                    EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::rook>(r, 9)), (grid_type::neighbor_indices_type {             { r-1, 9 },
                                                                                                                                  { r,   8 },  /*(r,c)*/  { r,   0 },
                                                                                                                                              { r+1, 9 }              }));
                }
            }

            { // Last row 
                // Bottom-left corner
                EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::queen>(4, 0)), (grid_type::neighbor_indices_type { { 3, 9 }, { 3, 0 }, { 3, 1 },
                                                                                                                               { 4, 9 }, /*(r,c)*/ { 4, 1 },
                                                                                                                               { 0, 9 }, { 0, 0 }, { 0, 1 } }));
                EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::rook>(4, 0)), (grid_type::neighbor_indices_type {           { 3, 0 },
                                                                                                                              { 4, 9 }, /*(r,c)*/ { 4, 1 },
                                                                                                                                        { 0, 0 }            }));
                // Inner cols
                for(std::size_t c=1; c<9; ++c)
                {
                    EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::queen>(4, c)), (grid_type::neighbor_indices_type { { 3, c-1 }, { 3, c }, { 3, c+1 },
                                                                                                                                   { 4, c-1 }, /*(r,c)*/ { 4, c+1 },
                                                                                                                                   { 0, c-1 }, { 0, c }, { 0, c+1 } }));
                    EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::rook>(4, c)), (grid_type::neighbor_indices_type {             { 3, c },
                                                                                                                                  { 4, c-1 }, /*(r,c)*/ { 4, c+1 },
                                                                                                                                              { 0, c }              }));
                }
                // Bottom-right corner
                EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::queen>(4, 9)), (grid_type::neighbor_indices_type { { 3, 8 }, { 3, 9 }, { 3, 0 },
                                                                                                                               { 4, 8 }, /*(r,c)*/ { 4, 0 },
                                                                                                                               { 0, 8 }, { 0, 9 }, { 0, 0 } }));
                EXPECT_EQ((looped_grid.neighbor_indices<fs::raster_connect::rook>(4, 9)), (grid_type::neighbor_indices_type {           { 3, 9 },
                                                                                                                              { 4, 8 }, /*(r,c)*/ { 4, 0 },
                                                                                                                                        { 0, 9 }            }));
            }
        }

        TEST_F(raster_grid, spacing)
        {
            EXPECT_EQ(fixed_grid.spacing(), (std::array<double, 2> {{1.3, 1.2}}));
            EXPECT_EQ(looped_grid.spacing(), (std::array<double, 2> {{1.4, 1.8}}));
        }

        TEST_F(raster_grid, size)
        {
            EXPECT_EQ(fixed_grid.size(), 50);
            EXPECT_EQ(looped_grid.size(), 50);
        }

    }
}