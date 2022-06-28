#include "test_raster_grid.hpp"


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace testing
    {

        TEST_F(bishop_raster_grid, n_neighbors_max)
        {
            EXPECT_EQ(grid_type::n_neighbors_max(), 4u);
        }

#define EXPECT_INDICES(ROW, COL, INDICES, RC_INDICES)                                              \
    EXPECT_EQ(bishop_fixed.neighbors_indices(ROW* shape[1] + COL),                                 \
              (xt::xtensor<std::size_t, 1> INDICES));                                              \
                                                                                                   \
    bishop_fixed.neighbors_indices(ROW* shape[1] + COL, neighbors_idx);                            \
    EXPECT_EQ(neighbors_idx, (xt::xtensor<std::size_t, 1> INDICES));                               \
                                                                                                   \
    EXPECT_EQ((bishop_fixed.neighbors_indices(ROW, COL)),                                          \
              (grid_type::neighbors_indices_raster_type RC_INDICES));

        TEST_F(bishop_raster_grid_fixed, neighbors_indices)
        {
            {  // First row
                // Top-left corner
                EXPECT_INDICES(0,
                               0,
                               ({ 11 }),
                               ({ /*(r,c)*/
                                  { 1, 1 } }));
                // Inner cols
                for (std::size_t c = 1; c < 9; ++c)
                {
                    EXPECT_INDICES(0,
                                   c,
                                   ({ shape[1] + c - 1, shape[1] + c + 1 }),
                                   ({ /*(r,c)*/
                                      { 1, c - 1 },
                                      { 1, c + 1 } }));
                }
                // Top-right corner
                EXPECT_INDICES(0,
                               9,
                               ({ 18 }),
                               ({
                                   /*(r,c)*/
                                   { 1, 8 },
                               }));
            }
            {  // Inners rows
                for (std::size_t r = 1; r < 4; ++r)
                {
                    // First col
                    EXPECT_INDICES(r,
                                   0,
                                   ({ (r - 1) * shape[1] + 1, (r + 1) * shape[1] + 1 }),
                                   ({ { r - 1, 1 },
                                      /*(r,c)*/
                                      { r + 1, 1 } }));
                    // Inners cols
                    for (std::size_t c = 1; c < 9; ++c)
                    {
                        EXPECT_INDICES(r,
                                       c,
                                       ({ (r - 1) * shape[1] + c - 1,
                                          (r - 1) * shape[1] + c + 1,
                                          (r + 1) * shape[1] + c - 1,
                                          (r + 1) * shape[1] + c + 1 }),
                                       ({ { r - 1, c - 1 },
                                          { r - 1, c + 1 },
                                          /*(r,c)*/
                                          { r + 1, c - 1 },
                                          { r + 1, c + 1 } }));
                    }
                    // last col
                    EXPECT_INDICES(r,
                                   shape[1] - 1,
                                   ({ (r - 1) * shape[1] + 8, (r + 1) * shape[1] + 8 }),
                                   ({
                                       { r - 1, 8 },
                                       /*(r,c)*/
                                       { r + 1, 8 },
                                   }));
                }
            }

            {  // Last row
                // Bottom-left corner
                EXPECT_INDICES(4,
                               0,
                               ({ 31 }),
                               ({ { 3, 1 },
                                  /*(r,c)*/ }));
                // Inner cols
                for (std::size_t c = 1; c < 9; ++c)
                {
                    EXPECT_INDICES(4,
                                   c,
                                   ({ 30 + c - 1, 30 + c + 1 }),
                                   ({ { 3, c - 1 }, { 3, c + 1 },
                                      /*(r,c)*/ }));
                }
                // Bottom-right corner
                EXPECT_INDICES(4,
                               9,
                               ({ 38 }),
                               ({ { 3, 8 }
                                  /*(r,c)*/ }));
            }
        }
#undef EXPECT_INDICES

#define EXPECT_INDICES(ROW, COL, INDICES, RC_INDICES)                                              \
    EXPECT_EQ(bishop_looped.neighbors_indices(ROW* shape[1] + COL),                                \
              (xt::xtensor<std::size_t, 1> INDICES));                                              \
                                                                                                   \
    bishop_looped.neighbors_indices(ROW* shape[1] + COL, neighbors_idx);                           \
    EXPECT_EQ(neighbors_idx, (xt::xtensor<std::size_t, 1> INDICES));                               \
                                                                                                   \
    EXPECT_EQ((bishop_looped.neighbors_indices(ROW, COL)),                                         \
              (grid_type::neighbors_indices_raster_type RC_INDICES));

        TEST_F(bishop_raster_grid_looped, neighbors_indices)
        {
            {  // First row
                // Top-left corner
                EXPECT_INDICES(0,
                               0,
                               ({ 49, 41, 19, 11 }),
                               ({ { 4, 9 },
                                  { 4, 1 },
                                  /*(r,c)*/
                                  { 1, 9 },
                                  { 1, 1 } }));
                // Inner cols
                for (std::size_t c = 1; c < 9; ++c)
                {
                    EXPECT_INDICES(0,
                                   c,
                                   ({ 40 + c - 1, 40 + c + 1, 10 + c - 1, 10 + c + 1 }),
                                   ({ { 4, c - 1 },
                                      { 4, c + 1 },
                                      /*(r,c)*/
                                      { 1, c - 1 },
                                      { 1, c + 1 } }));
                }
                // Top-right corner
                EXPECT_INDICES(0,
                               9,
                               ({ 48, 40, 18, 10 }),
                               ({ { 4, 8 },
                                  { 4, 0 },
                                  /*(r,c)*/
                                  { 1, 8 },
                                  { 1, 0 } }));
            }
            {  // Inners rows
                for (std::size_t r = 1; r < 4; ++r)
                {
                    // First col
                    EXPECT_INDICES(r,
                                   0,
                                   ({ (r - 1) * shape[1] + 9,
                                      (r - 1) * shape[1] + 1,
                                      (r + 1) * shape[1] + 9,
                                      (r + 1) * shape[1] + 1 }),
                                   ({ { r - 1, 9 },
                                      { r - 1, 1 },
                                      /*(r,c)*/
                                      { r + 1, 9 },
                                      { r + 1, 1 } }));
                    // Inners cols
                    for (std::size_t c = 1; c < 9; ++c)
                    {
                        EXPECT_INDICES(r,
                                       c,
                                       ({ (r - 1) * shape[1] + c - 1,
                                          (r - 1) * shape[1] + c + 1,
                                          (r + 1) * shape[1] + c - 1,
                                          (r + 1) * shape[1] + c + 1 }),
                                       ({ { r - 1, c - 1 },
                                          { r - 1, c + 1 },
                                          /*(r,c)*/
                                          { r + 1, c - 1 },
                                          { r + 1, c + 1 } }));
                    }
                    // last col
                    EXPECT_INDICES(r,
                                   shape[1] - 1,
                                   ({ (r - 1) * shape[1] + 8,
                                      (r - 1) * shape[1],
                                      (r + 1) * shape[1] + 8,
                                      (r + 1) * shape[1] }),
                                   ({ { r - 1, 8 },
                                      { r - 1, 0 },
                                      /*(r,c)*/
                                      { r + 1, 8 },
                                      { r + 1, 0 } }));
                }
            }
            {  // Last row
                // Bottom-left corner
                EXPECT_INDICES(4,
                               0,
                               ({ 39, 31, 9, 1 }),
                               ({ { 3, 9 },
                                  { 3, 1 },
                                  /*(r,c)*/
                                  { 0, 9 },
                                  { 0, 1 } }));
                // Inner cols
                for (std::size_t c = 1; c < 9; ++c)
                {
                    EXPECT_INDICES(4,
                                   c,
                                   ({ 3 * shape[1] + c - 1, 3 * shape[1] + c + 1, c - 1, c + 1 }),
                                   ({ { 3, c - 1 },
                                      { 3, c + 1 },
                                      /*(r,c)*/
                                      { 0, c - 1 },
                                      { 0, c + 1 } }));
                }
                // Bottom-right corner
                EXPECT_INDICES(4,
                               9,
                               ({ 38, 30, 8, 0 }),
                               ({ { 3, 8 },
                                  { 3, 0 },
                                  /*(r,c)*/
                                  { 0, 8 },
                                  { 0, 0 } }));
            }
        }
#undef EXPECT_INDICES


        TEST_F(bishop_raster_grid_fixed, neighbors_distances)
        {
            using grid_data_type = grid_type::grid_data_type;
            using neighbors_distances_type = grid_type::neighbors_distances_type;

            grid_data_type dia = std::sqrt(1.2 * 1.2 + 1.3 * 1.3);

            // Top-left corner
            EXPECT_EQ(bishop_fixed.neighbors_distances(0), (neighbors_distances_type{ dia }));

            // Bottom-right corner
            EXPECT_EQ(bishop_fixed.neighbors_distances(49), (neighbors_distances_type{ dia }));

            // Inners rows/cols
            for (std::size_t r = 1; r < 4; ++r)
            {
                for (std::size_t c = 1; c < 9; ++c)
                {
                    EXPECT_EQ(bishop_fixed.neighbors_distances(r * shape[1] + c),
                              (neighbors_distances_type{ dia, dia, dia, dia }));
                }
            }
        }

        TEST_F(bishop_raster_grid_looped, neighbors_distances)
        {
            using grid_data_type = grid_type::grid_data_type;
            using neighbors_distances_type = grid_type::neighbors_distances_type;

            grid_data_type dia = std::sqrt(1.4 * 1.4 + 1.8 * 1.8);

            for (std::size_t r = 0; r < 5; ++r)
            {
                for (std::size_t c = 0; c < 10; ++c)
                {
                    EXPECT_EQ(bishop_looped.neighbors_distances(r * shape[1] + c),
                              (neighbors_distances_type{ dia, dia, dia, dia }));
                }
            }
        }

        TEST_F(bishop_raster_grid_fixed, neighbors)
        {
            using neighbors_type = std::vector<fs::neighbor>;

            double d1 = std::sqrt(1.3 * 1.3 + 1.2 * 1.2);
            double d2 = std::sqrt(1.4 * 1.4 + 1.8 * 1.8);

            {  // First row
                // Top-left corner
                EXPECT_EQ(bishop_fixed.neighbors(0),
                          (neighbors_type{ /* Node */
                                           { 11, d1, co } }));

                EXPECT_EQ(bishop_looped.neighbors(0),
                          (neighbors_type{ { 49, d2, lb },
                                           { 41, d2, lb },
                                           /* Node */
                                           { 19, d2, lb },
                                           { 11, d2, co } }));

                // Inner cols
                for (std::size_t c = 2; c < 8; ++c)
                {
                    EXPECT_EQ(bishop_fixed.neighbors(c),
                              (neighbors_type{ /* Node */
                                               { c + 9, d1, co },
                                               { c + 11, d1, co } }));

                    EXPECT_EQ(bishop_looped.neighbors(c),
                              (neighbors_type{ { c + 39, d2, lb },
                                               { c + 41, d2, lb },
                                               /* Node */
                                               { c + 9, d2, co },
                                               { c + 11, d2, co } }));
                }

                EXPECT_EQ(bishop_looped.neighbors(22),
                          (neighbors_type{ { 11, d2, co },
                                           { 13, d2, co },
                                           /* Node */
                                           { 31, d2, co },
                                           { 33, d2, co } }));
            }
        }

        TEST_F(bishop_raster_grid, neighbors_count)
        {
            // Top-left corner nodes
            EXPECT_EQ(bishop_fixed.neighbors_count(static_cast<std::size_t>(0)), 1);
            EXPECT_EQ(bishop_hlooped.neighbors_count(static_cast<std::size_t>(0)), 2);
            EXPECT_EQ(bishop_vlooped.neighbors_count(static_cast<std::size_t>(0)), 2);
            EXPECT_EQ(bishop_looped.neighbors_count(static_cast<std::size_t>(0)), 4);

            // Top-right corner nodes
            EXPECT_EQ(bishop_fixed.neighbors_count(shape[1] - 1), 1);
            EXPECT_EQ(bishop_hlooped.neighbors_count(shape[1] - 1), 2);
            EXPECT_EQ(bishop_vlooped.neighbors_count(shape[1] - 1), 2);
            EXPECT_EQ(bishop_looped.neighbors_count(shape[1] - 1), 4);

            // Bottom-left corner nodes
            EXPECT_EQ(bishop_fixed.neighbors_count((shape[0] - 1) * shape[1]), 1);
            EXPECT_EQ(bishop_hlooped.neighbors_count((shape[0] - 1) * shape[1]), 2);
            EXPECT_EQ(bishop_vlooped.neighbors_count((shape[0] - 1) * shape[1]), 2);
            EXPECT_EQ(bishop_looped.neighbors_count((shape[0] - 1) * shape[1]), 4);

            // Bottom-right corner nodes
            EXPECT_EQ(bishop_fixed.neighbors_count(shape[0] * shape[1] - 1), 1);
            EXPECT_EQ(bishop_hlooped.neighbors_count(shape[0] * shape[1] - 1), 2);
            EXPECT_EQ(bishop_vlooped.neighbors_count(shape[0] * shape[1] - 1), 2);
            EXPECT_EQ(bishop_looped.neighbors_count(shape[0] * shape[1] - 1), 4);

            for (std::size_t c = 1; c < shape[1] - 1; ++c)
            {
                // Top edge nodes (without corners)
                EXPECT_EQ(bishop_fixed.neighbors_count(c), 2);
                EXPECT_EQ(bishop_hlooped.neighbors_count(c), 2);
                EXPECT_EQ(bishop_vlooped.neighbors_count(c), 4);
                EXPECT_EQ(bishop_looped.neighbors_count(c), 4);

                // Bottom edge nodes (without corners)
                EXPECT_EQ(bishop_fixed.neighbors_count((shape[0] - 1) * shape[1] + c), 2);
                EXPECT_EQ(bishop_hlooped.neighbors_count((shape[0] - 1) * shape[1] + c), 2);
                EXPECT_EQ(bishop_vlooped.neighbors_count((shape[0] - 1) * shape[1] + c), 4);
                EXPECT_EQ(bishop_looped.neighbors_count((shape[0] - 1) * shape[1] + c), 4);
            }

            for (std::size_t r = 1; r < shape[0] - 1; ++r)
            {
                // Left edge nodes (without corners)
                EXPECT_EQ(bishop_fixed.neighbors_count(r * shape[1]), 2);
                EXPECT_EQ(bishop_hlooped.neighbors_count(r * shape[1]), 4);
                EXPECT_EQ(bishop_vlooped.neighbors_count(r * shape[1]), 2);
                EXPECT_EQ(bishop_looped.neighbors_count(r * shape[1]), 4);

                // Right edge nodes (without corners)
                EXPECT_EQ(bishop_fixed.neighbors_count((r + 1) * shape[1] - 1), 2);
                EXPECT_EQ(bishop_hlooped.neighbors_count((r + 1) * shape[1] - 1), 4);
                EXPECT_EQ(bishop_vlooped.neighbors_count((r + 1) * shape[1] - 1), 2);
                EXPECT_EQ(bishop_looped.neighbors_count((r + 1) * shape[1] - 1), 4);
            }

            for (std::size_t r = 1; r < shape[0] - 1; ++r)
            {
                for (std::size_t c = 1; c < shape[1] - 1; ++c)
                {
                    // Inner nodes
                    EXPECT_EQ(bishop_fixed.neighbors_count(r * shape[1] + c), 4);
                    EXPECT_EQ(bishop_hlooped.neighbors_count(r * shape[1] + c), 4);
                    EXPECT_EQ(bishop_vlooped.neighbors_count(r * shape[1] + c), 4);
                    EXPECT_EQ(bishop_looped.neighbors_count(r * shape[1] + c), 4);
                }
            }
        }
    }
}
