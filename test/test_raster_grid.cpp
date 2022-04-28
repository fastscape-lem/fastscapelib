#include "test_raster_grid.hpp"


namespace fs = fastscapelib;


namespace fastscapelib
{
    namespace testing
    {

        TEST_F(raster_boundary_status, ctor)
        {
            EXPECT_EQ(fixed_value_status.left, fb);
            EXPECT_EQ(fixed_value_status.right, fb);
            EXPECT_EQ(fixed_value_status.bottom, fb);
            EXPECT_EQ(fixed_value_status.top, fb);

            EXPECT_EQ(hlooped_status.left, lb);
            EXPECT_EQ(hlooped_status.right, lb);
            EXPECT_EQ(hlooped_status.bottom, fb);
            EXPECT_EQ(hlooped_status.top, fb);

            EXPECT_EQ(vlooped_status.left, fb);
            EXPECT_EQ(vlooped_status.right, fb);
            EXPECT_EQ(vlooped_status.bottom, lb);
            EXPECT_EQ(vlooped_status.top, lb);

            EXPECT_EQ(hvlooped_status.left, lb);
            EXPECT_EQ(hvlooped_status.right, lb);
            EXPECT_EQ(hvlooped_status.bottom, lb);
            EXPECT_EQ(hvlooped_status.top, lb);

            EXPECT_THROW(fs::raster_boundary_status{ hill_formed_loop }, std::invalid_argument);
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

        TEST_F(raster_grid, ctor)
        {
            node_status lg = node_status::fixed_gradient_boundary;

            std::array<size_type, 2> shape{ { 3, 3 } };
            std::vector<fs::raster_node> nodes_vector1{ fs::raster_node(
                { 1, 1, fs::node_status::fixed_gradient_boundary }) };
            auto g1 = grid_type(shape, { 1.4, 1.8 }, hloop, nodes_vector1);

            auto expected_status
                = grid_type::node_status_type{ { fb, fb, fb }, { lb, lg, lb }, { fb, fb, fb } };
            ASSERT_EQ(g1.status_at_nodes(), expected_status);

            std::vector<fs::raster_node> nodes_vector2{ fs::raster_node({ 15, 15, co }) };
            ASSERT_THROW(grid_type(shape, { 1.4, 1.8 }, fb, nodes_vector2), std::out_of_range);

            std::vector<fs::raster_node> nodes_vector3{ fs::raster_node({ 0, 0, co }) };
            ASSERT_THROW(grid_type(shape, { 1.4, 1.8 }, hvloop, nodes_vector3),
                         std::invalid_argument);
        }

        TEST_F(raster_grid, static_expr)
        {
            EXPECT_EQ(fs::raster_grid::is_structured(), true);
            EXPECT_EQ(fs::raster_grid::is_uniform(), true);
            EXPECT_EQ(fs::raster_grid::max_neighbors(), 8u);
            EXPECT_EQ(fs::raster_grid::xt_ndims(), 2);
        }

        TEST_F(raster_grid, spacing)
        {
            EXPECT_TRUE(xt::all(xt::equal(fixed_grid.spacing(), spacing_type({ 1.3, 1.2 }))));
            EXPECT_TRUE(xt::all(xt::equal(looped_grid.spacing(), spacing_type({ 1.4, 1.8 }))));
        }

        TEST_F(raster_grid, size)
        {
            EXPECT_EQ(fixed_grid.size(), 50u);
            EXPECT_EQ(looped_grid.size(), 50u);
        }

        TEST_F(raster_grid, length)
        {
            EXPECT_TRUE(xt::all(xt::isclose(fixed_grid.length(), length_type({ 5.2, 10.8 }))));
            EXPECT_TRUE(xt::all(xt::equal(looped_grid.length(), length_type({ 5.6, 16.2 }))));
        }

        TEST_F(raster_grid, node_area)
        {
            for (auto n : fixed_grid.nodes_indices())
            {
                EXPECT_EQ(fixed_grid.node_area(n), 1.56);
                EXPECT_EQ(looped_grid.node_area(n), 2.52);
            }
        }

        TEST_F(raster_grid, from_length)
        {
            auto grid_from_length = grid_type::from_length(
                shape_type({ { 151, 101 } }), length_type({ 1500., 2000. }), fixed_value_status);
            EXPECT_TRUE(
                xt::all(xt::equal(grid_from_length.length(), length_type({ 1500., 2000. }))));
            EXPECT_EQ(grid_from_length.size(), 15251u);
            EXPECT_TRUE(xt::all(xt::equal(grid_from_length.spacing(), spacing_type({ 10., 20. }))));
        }

        TEST_F(raster_grid, clone)
        {
            using queen_type = fs::raster_grid_xt<fs::xt_selector, fs::raster_connect::queen>;
            using neighbors_type = std::vector<fs::neighbor>;

            double d1 = std::sqrt(1.3 * 1.3 + 1.2 * 1.2);

            auto queen_fixed = queen_type(shape, { 1.3, 1.2 }, fixed_value_status);
            EXPECT_EQ(queen_fixed.neighbors(9),        // Top-right corner
                      (neighbors_type{ { 8, 1.2, fb }, /* Node */
                                       { 18, d1, co },
                                       { 19, 1.3, fb } }));

            auto rook_like = fs::raster_grid::clone<raster_connect::rook>(queen_fixed);
            EXPECT_EQ(rook_like.neighbors(9),          // Top-right corner
                      (neighbors_type{ { 8, 1.2, fb }, /* Node */
                                       { 19, 1.3, fb } }));

            auto bishop_like = fs::raster_grid::clone<raster_connect::bishop>(queen_fixed);
            EXPECT_EQ(bishop_like.neighbors(9),  // Top-right corner
                      (neighbors_type{           /* Node */
                                       { 18, d1, co } }));
        }

        TEST_F(raster_grid, node_code)
        {
            // Top-left corner nodes
            EXPECT_EQ(fixed_grid.node_code(0), 0);
            EXPECT_EQ(looped_grid.node_code(0), 0);

            // Top-right corner nodes
            EXPECT_EQ(fixed_grid.node_code(shape[1] - 1), 2);
            EXPECT_EQ(looped_grid.node_code(shape[1] - 1), 2);

            // Bottom-left corner nodes
            EXPECT_EQ(fixed_grid.node_code((shape[0] - 1) * shape[1]), 6);
            EXPECT_EQ(looped_grid.node_code((shape[0] - 1) * shape[1]), 6);

            // Bottom-right corner nodes
            EXPECT_EQ(fixed_grid.node_code(shape[0] * shape[1] - 1), 8);
            EXPECT_EQ(looped_grid.node_code(shape[0] * shape[1] - 1), 8);

            for (std::size_t c = 1; c < shape[1] - 1; ++c)
            {
                // Top edge nodes (without corners)
                EXPECT_EQ(fixed_grid.node_code(c), 1);
                EXPECT_EQ(looped_grid.node_code(c), 1);

                // Bottom edge nodes (without corners)
                EXPECT_EQ(fixed_grid.node_code((shape[0] - 1) * shape[1] + c), 7);
                EXPECT_EQ(looped_grid.node_code((shape[0] - 1) * shape[1] + c), 7);
            }

            for (std::size_t r = 1; r < shape[0] - 1; ++r)
            {
                // Left edge nodes (without corners)
                EXPECT_EQ(fixed_grid.node_code(r * shape[1]), 3);
                EXPECT_EQ(looped_grid.node_code(r * shape[1]), 3);

                // Right edge nodes (without corners)
                EXPECT_EQ(fixed_grid.node_code((r + 1) * shape[1] - 1), 5);
                EXPECT_EQ(looped_grid.node_code((r + 1) * shape[1] - 1), 5);
            }

            for (std::size_t r = 1; r < shape[0] - 1; ++r)
            {
                for (std::size_t c = 1; c < shape[1] - 1; ++c)
                {
                    // Inner nodes
                    EXPECT_EQ(fixed_grid.node_code(r * 10 + c), 4);
                    EXPECT_EQ(looped_grid.node_code(r * 10 + c), 4);
                }
            }
        }
    }
}
