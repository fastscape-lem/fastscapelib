#include "test_raster_grid.hpp"


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace testing
    {

        TEST_F(queen_raster_grid, n_neighbors_max)
        {
            EXPECT_EQ(grid_type::n_neighbors_max(), 8u);
        }

#define EXPECT_INDICES(ROW, COL, INDICES, RC_INDICES)                                              \
    EXPECT_EQ(queen_fixed.neighbors_indices(ROW* shape[1] + COL),                                  \
              (xt::xtensor<std::size_t, 1> INDICES));                                              \
                                                                                                   \
    queen_fixed.neighbors_indices(ROW* shape[1] + COL, neighbors_idx);                             \
    EXPECT_EQ(neighbors_idx, (xt::xtensor<std::size_t, 1> INDICES));                               \
                                                                                                   \
    EXPECT_EQ((queen_fixed.neighbors_indices(ROW, COL)),                                           \
              (grid_type::neighbors_indices_raster_type RC_INDICES));

        TEST_F(queen_raster_grid_fixed, neighbors_indices)
        {
            {  // First row
                // Top-left corner
                EXPECT_INDICES(0, 0, ({ 1, 10, 11 }), ({ /*(r,c)*/ { 0, 1 }, { 1, 0 }, { 1, 1 } }));
                // Inner cols
                for (std::size_t c = 1; c < 9; ++c)
                {
                    EXPECT_INDICES(
                        0,
                        c,
                        ({ c - 1, c + 1, shape[1] + c - 1, shape[1] + c, shape[1] + c + 1 }),
                        ({ { 0, c - 1 },
                           /*(r,c)*/ { 0, c + 1 },
                           { 1, c - 1 },
                           { 1, c },
                           { 1, c + 1 } }));
                }
                // Top-right corner
                EXPECT_INDICES(0,
                               9,
                               ({ 8, 18, 19 }),
                               ({ { 0, 8 }, /*(r,c)*/
                                  { 1, 8 },
                                  { 1, 9 } }));
            }
            {  // Inners rows
                for (std::size_t r = 1; r < 4; ++r)
                {
                    // First col
                    EXPECT_INDICES(r,
                                   0,
                                   ({ (r - 1) * shape[1],
                                      (r - 1) * shape[1] + 1,
                                      r * shape[1] + 1,
                                      (r + 1) * shape[1],
                                      (r + 1) * shape[1] + 1 }),
                                   ({ { r - 1, 0 },
                                      { r - 1, 1 },
                                      /*(r,c)*/ { r, 1 },
                                      { r + 1, 0 },
                                      { r + 1, 1 } }));
                    // Inners cols
                    for (std::size_t c = 1; c < 9; ++c)
                    {
                        EXPECT_INDICES(r,
                                       c,
                                       ({ (r - 1) * shape[1] + c - 1,
                                          (r - 1) * shape[1] + c,
                                          (r - 1) * shape[1] + c + 1,
                                          r * shape[1] + c - 1,
                                          r * shape[1] + c + 1,
                                          (r + 1) * shape[1] + c - 1,
                                          (r + 1) * shape[1] + c,
                                          (r + 1) * shape[1] + c + 1 }),
                                       ({ { r - 1, c - 1 },
                                          { r - 1, c },
                                          { r - 1, c + 1 },
                                          { r, c - 1 },
                                          /*(r,c)*/ { r, c + 1 },
                                          { r + 1, c - 1 },
                                          { r + 1, c },
                                          { r + 1, c + 1 } }));
                    }
                    // last col
                    EXPECT_INDICES(r,
                                   shape[1] - 1,
                                   ({ (r - 1) * shape[1] + 8,
                                      (r - 1) * shape[1] + 9,
                                      r * shape[1] + 8,
                                      (r + 1) * shape[1] + 8,
                                      (r + 1) * shape[1] + 9 }),
                                   ({ { r - 1, 8 },
                                      { r - 1, 9 },
                                      { r, 8 }, /*(r,c)*/
                                      { r + 1, 8 },
                                      { r + 1, 9 } }));
                }
            }
            {  // Last row
                // Bottom-left corner
                EXPECT_INDICES(4,
                               0,
                               ({ 30, 31, 41 }),
                               ({ { 3, 0 },
                                  { 3, 1 },
                                  /*(r,c)*/ { 4, 1 } }));
                // Inner cols
                for (std::size_t c = 1; c < 9; ++c)
                {
                    EXPECT_INDICES(4,
                                   c,
                                   ({ 30 + c - 1, 30 + c, 30 + c + 1, 40 + c - 1, 40 + c + 1 }),
                                   ({ { 3, c - 1 },
                                      { 3, c },
                                      { 3, c + 1 },
                                      { 4, c - 1 },
                                      /*(r,c)*/ { 4, c + 1 } }));
                }
                // Bottom-right corner
                EXPECT_INDICES(
                    4, 9, ({ 38, 39, 48 }), ({ { 3, 8 }, { 3, 9 }, { 4, 8 } /*(r,c)*/ }));
            }
        }
#undef EXPECT_INDICES

#define EXPECT_INDICES(ROW, COL, INDICES, RC_INDICES)                                              \
    EXPECT_EQ(queen_looped.neighbors_indices(ROW* shape[1] + COL),                                 \
              (xt::xtensor<std::size_t, 1> INDICES));                                              \
                                                                                                   \
    queen_looped.neighbors_indices(ROW* shape[1] + COL, neighbors_idx);                            \
    EXPECT_EQ(neighbors_idx, (xt::xtensor<std::size_t, 1> INDICES));                               \
                                                                                                   \
    EXPECT_EQ((queen_looped.neighbors_indices(ROW, COL)),                                          \
              (grid_type::neighbors_indices_raster_type RC_INDICES));

        TEST_F(queen_raster_grid_looped, neighbors_indices)
        {
            {  // First row
                // Top-left corner
                EXPECT_INDICES(0,
                               0,
                               ({ 49, 40, 41, 9, 1, 19, 10, 11 }),
                               ({ { 4, 9 },
                                  { 4, 0 },
                                  { 4, 1 },
                                  { 0, 9 },
                                  /*(r,c)*/ { 0, 1 },
                                  { 1, 9 },
                                  { 1, 0 },
                                  { 1, 1 } }));
                // Inner cols
                for (std::size_t c = 1; c < 9; ++c)
                {
                    EXPECT_INDICES(0,
                                   c,
                                   ({ 40 + c - 1,
                                      40 + c,
                                      40 + c + 1,
                                      c - 1,
                                      c + 1,
                                      10 + c - 1,
                                      10 + c,
                                      10 + c + 1 }),
                                   ({ { 4, c - 1 },
                                      { 4, c },
                                      { 4, c + 1 },
                                      { 0, c - 1 },
                                      /*(r,c)*/ { 0, c + 1 },
                                      { 1, c - 1 },
                                      { 1, c },
                                      { 1, c + 1 } }));
                }
                // Top-right corner
                EXPECT_INDICES(0,
                               9,
                               ({ 48, 49, 40, 8, 0, 18, 19, 10 }),
                               ({ { 4, 8 },
                                  { 4, 9 },
                                  { 4, 0 },
                                  { 0, 8 },
                                  /*(r,c)*/ { 0, 0 },
                                  { 1, 8 },
                                  { 1, 9 },
                                  { 1, 0 } }));
            }
            {  // Inners rows
                for (std::size_t r = 1; r < 4; ++r)
                {
                    // First col
                    EXPECT_INDICES(r,
                                   0,
                                   ({ (r - 1) * shape[1] + 9,
                                      (r - 1) * shape[1],
                                      (r - 1) * shape[1] + 1,
                                      r * shape[1] + 9,
                                      r * shape[1] + 1,
                                      (r + 1) * shape[1] + 9,
                                      (r + 1) * shape[1],
                                      (r + 1) * shape[1] + 1 }),
                                   ({ { r - 1, 9 },
                                      { r - 1, 0 },
                                      { r - 1, 1 },
                                      { r, 9 },
                                      /*(r,c)*/ { r, 1 },
                                      { r + 1, 9 },
                                      { r + 1, 0 },
                                      { r + 1, 1 } }));
                    // Inners cols
                    for (std::size_t c = 1; c < 9; ++c)
                    {
                        EXPECT_INDICES(r,
                                       c,
                                       ({ (r - 1) * shape[1] + c - 1,
                                          (r - 1) * shape[1] + c,
                                          (r - 1) * shape[1] + c + 1,
                                          r * shape[1] + c - 1,
                                          r * shape[1] + c + 1,
                                          (r + 1) * shape[1] + c - 1,
                                          (r + 1) * shape[1] + c,
                                          (r + 1) * shape[1] + c + 1 }),
                                       ({ { r - 1, c - 1 },
                                          { r - 1, c },
                                          { r - 1, c + 1 },
                                          { r, c - 1 },
                                          /*(r,c)*/ { r, c + 1 },
                                          { r + 1, c - 1 },
                                          { r + 1, c },
                                          { r + 1, c + 1 } }));
                    }
                    // last col
                    EXPECT_INDICES(r,
                                   shape[1] - 1,
                                   ({ (r - 1) * shape[1] + 8,
                                      (r - 1) * shape[1] + 9,
                                      (r - 1) * shape[1],
                                      r * shape[1] + 8,
                                      r * shape[1],
                                      (r + 1) * shape[1] + 8,
                                      (r + 1) * shape[1] + 9,
                                      (r + 1) * shape[1] }),
                                   ({ { r - 1, 8 },
                                      { r - 1, 9 },
                                      { r - 1, 0 },
                                      { r, 8 },
                                      /*(r,c)*/ { r, 0 },
                                      { r + 1, 8 },
                                      { r + 1, 9 },
                                      { r + 1, 0 } }));
                }
            }
            {  // Last row
                // Bottom-left corner
                EXPECT_INDICES(4,
                               0,
                               ({ 39, 30, 31, 49, 41, 9, 0, 1 }),
                               ({ { 3, 9 },
                                  { 3, 0 },
                                  { 3, 1 },
                                  { 4, 9 },
                                  /*(r,c)*/ { 4, 1 },
                                  { 0, 9 },
                                  { 0, 0 },
                                  { 0, 1 } }));
                // Inner cols
                for (std::size_t c = 1; c < 9; ++c)
                {
                    EXPECT_INDICES(4,
                                   c,
                                   ({ 3 * shape[1] + c - 1,
                                      3 * shape[1] + c,
                                      3 * shape[1] + c + 1,
                                      4 * shape[1] + c - 1,
                                      4 * shape[1] + c + 1,
                                      c - 1,
                                      c,
                                      c + 1 }),
                                   ({ { 3, c - 1 },
                                      { 3, c },
                                      { 3, c + 1 },
                                      { 4, c - 1 },
                                      /*(r,c)*/ { 4, c + 1 },
                                      { 0, c - 1 },
                                      { 0, c },
                                      { 0, c + 1 } }));
                }
                // Bottom-right corner
                EXPECT_INDICES(4,
                               9,
                               ({ 38, 39, 30, 48, 40, 8, 9, 0 }),
                               ({ { 3, 8 },
                                  { 3, 9 },
                                  { 3, 0 },
                                  { 4, 8 },
                                  /*(r,c)*/ { 4, 0 },
                                  { 0, 8 },
                                  { 0, 9 },
                                  { 0, 0 } }));
            }
        }
#undef EXPECT_INDICES

        TEST_F(queen_raster_grid_fixed, neighbors_distances)
        {
            using grid_data_type = grid_type::grid_data_type;
            using neighbors_distances_type = grid_type::neighbors_distances_type;

            grid_data_type dia = std::sqrt(1.2 * 1.2 + 1.3 * 1.3);

            // Top-left corner
            EXPECT_EQ(queen_fixed.neighbors_distances(0),
                      (neighbors_distances_type{ 1.2, 1.3, dia }));

            // Bottom-right corner
            EXPECT_EQ(queen_fixed.neighbors_distances(49),
                      (neighbors_distances_type{ dia, 1.3, 1.2 }));

            // Inners rows/cols
            for (std::size_t r = 1; r < 4; ++r)
            {
                for (std::size_t c = 1; c < 9; ++c)
                {
                    EXPECT_EQ(queen_fixed.neighbors_distances(r * shape[1] + c),
                              (neighbors_distances_type{ dia, 1.3, dia, 1.2, 1.2, dia, 1.3, dia }));
                }
            }
        }

        TEST_F(queen_raster_grid_looped, neighbors_distances)
        {
            using grid_data_type = grid_type::grid_data_type;
            using neighbors_distances_type = grid_type::neighbors_distances_type;

            grid_data_type dia = std::sqrt(1.4 * 1.4 + 1.8 * 1.8);

            for (std::size_t r = 0; r < 5; ++r)
            {
                for (std::size_t c = 0; c < 10; ++c)
                {
                    EXPECT_EQ(queen_looped.neighbors_distances(r * shape[1] + c),
                              (neighbors_distances_type{ dia, 1.4, dia, 1.8, 1.8, dia, 1.4, dia }));
                }
            }
        }

        TEST_F(queen_raster_grid_fixed, neighbors)
        {  // Do not test extensively, indices and distances are already tested
            using neighbors_type = std::vector<fs::neighbor>;

            double d1 = std::sqrt(1.3 * 1.3 + 1.2 * 1.2);
            double d2 = std::sqrt(1.4 * 1.4 + 1.8 * 1.8);

            {  // First row
                // Top-left corner
                EXPECT_EQ(
                    queen_fixed.neighbors(0),
                    (neighbors_type{ /* Node */ { 1, 1.2, fb }, { 10, 1.3, fb }, { 11, d1, co } }));

                EXPECT_EQ(queen_looped.neighbors(0),
                          (neighbors_type{ { 49, d2, lb },
                                           { 40, 1.4, lb },
                                           { 41, d2, lb },
                                           { 9, 1.8, lb },
                                           /* Node */ { 1, 1.8, lb },
                                           { 19, d2, lb },
                                           { 10, 1.4, lb },
                                           { 11, d2, co } }));

                // Inner cols
                for (std::size_t c = 2; c < 8; ++c)
                {
                    EXPECT_EQ(queen_fixed.neighbors(c),
                              (neighbors_type{ { c - 1, 1.2, fb },
                                               /* Node */ { c + 1, 1.2, fb },
                                               { c + 9, d1, co },
                                               { c + 10, 1.3, co },
                                               { c + 11, d1, co } }));

                    EXPECT_EQ(queen_looped.neighbors(c),
                              (neighbors_type{ { c + 39, d2, lb },
                                               { c + 40, 1.4, lb },
                                               { c + 41, d2, lb },
                                               { c - 1, 1.8, lb },
                                               /* Node */ { c + 1, 1.8, lb },
                                               { c + 9, d2, co },
                                               { c + 10, 1.4, co },
                                               { c + 11, d2, co } }));
                }

                EXPECT_EQ(queen_looped.neighbors(22),
                          (neighbors_type{ { 11, d2, co },
                                           { 12, 1.4, co },
                                           { 13, d2, co },
                                           { 21, 1.8, co },
                                           /* Node */ { 23, 1.8, co },
                                           { 31, d2, co },
                                           { 32, 1.4, co },
                                           { 33, d2, co } }));
            }
        }

        TEST_F(queen_raster_grid_fixed, raster_neighbors)
        {
            using neighbors_type = std::vector<fs::raster_neighbor>;

            double d1 = std::sqrt(1.3 * 1.3 + 1.2 * 1.2);

            {
                EXPECT_EQ(queen_fixed.neighbors(0, 0),
                          (neighbors_type{ /* Node */ { 1, 0, 1, 1.2, fb },
                                           { 10, 1, 0, 1.3, fb },
                                           { 11, 1, 1, d1, co } }));

                EXPECT_EQ(queen_fixed.neighbors(1, 1),
                          (neighbors_type{ { 0, 0, 0, d1, fb },
                                           { 1, 0, 1, 1.3, fb },
                                           { 2, 0, 2, d1, fb },
                                           { 10, 1, 0, 1.2, fb },
                                           /* Node */ { 12, 1, 2, 1.2, co },
                                           { 20, 2, 0, d1, fb },
                                           { 21, 2, 1, 1.3, co },
                                           { 22, 2, 2, d1, co } }));
            }
        }

        TEST_F(queen_raster_grid, neighbors_count)
        {
            // Top-left corner nodes
            EXPECT_EQ(queen_fixed.neighbors_count(static_cast<std::size_t>(0)), 3);
            EXPECT_EQ(queen_hlooped.neighbors_count(static_cast<std::size_t>(0)), 5);
            EXPECT_EQ(queen_vlooped.neighbors_count(static_cast<std::size_t>(0)), 5);
            EXPECT_EQ(queen_looped.neighbors_count(static_cast<std::size_t>(0)), 8);

            // Top-right corner nodes
            EXPECT_EQ(queen_fixed.neighbors_count(shape[1] - 1), 3);
            EXPECT_EQ(queen_hlooped.neighbors_count(shape[1] - 1), 5);
            EXPECT_EQ(queen_vlooped.neighbors_count(shape[1] - 1), 5);
            EXPECT_EQ(queen_looped.neighbors_count(shape[1] - 1), 8);

            // Bottom-left corner nodes
            EXPECT_EQ(queen_fixed.neighbors_count((shape[0] - 1) * shape[1]), 3);
            EXPECT_EQ(queen_hlooped.neighbors_count((shape[0] - 1) * shape[1]), 5);
            EXPECT_EQ(queen_vlooped.neighbors_count((shape[0] - 1) * shape[1]), 5);
            EXPECT_EQ(queen_looped.neighbors_count((shape[0] - 1) * shape[1]), 8);

            // Bottom-right corner nodes
            EXPECT_EQ(queen_fixed.neighbors_count(shape[0] * shape[1] - 1), 3);
            EXPECT_EQ(queen_hlooped.neighbors_count(shape[0] * shape[1] - 1), 5);
            EXPECT_EQ(queen_vlooped.neighbors_count(shape[0] * shape[1] - 1), 5);
            EXPECT_EQ(queen_looped.neighbors_count(shape[0] * shape[1] - 1), 8);

            for (std::size_t c = 1; c < shape[1] - 1; ++c)
            {
                // Top edge nodes (without corners)
                EXPECT_EQ(queen_fixed.neighbors_count(c), 5);
                EXPECT_EQ(queen_hlooped.neighbors_count(c), 5);
                EXPECT_EQ(queen_vlooped.neighbors_count(c), 8);
                EXPECT_EQ(queen_looped.neighbors_count(c), 8);

                // Bottom edge nodes (without corners)
                EXPECT_EQ(queen_fixed.neighbors_count((shape[0] - 1) * shape[1] + c), 5);
                EXPECT_EQ(queen_hlooped.neighbors_count((shape[0] - 1) * shape[1] + c), 5);
                EXPECT_EQ(queen_vlooped.neighbors_count((shape[0] - 1) * shape[1] + c), 8);
                EXPECT_EQ(queen_looped.neighbors_count((shape[0] - 1) * shape[1] + c), 8);
            }

            for (std::size_t r = 1; r < shape[0] - 1; ++r)
            {
                // Left edge nodes (without corners)
                EXPECT_EQ(queen_fixed.neighbors_count(r * shape[1]), 5);
                EXPECT_EQ(queen_hlooped.neighbors_count(r * shape[1]), 8);
                EXPECT_EQ(queen_vlooped.neighbors_count(r * shape[1]), 5);
                EXPECT_EQ(queen_looped.neighbors_count(r * shape[1]), 8);

                // Right edge nodes (without corners)
                EXPECT_EQ(queen_fixed.neighbors_count((r + 1) * shape[1] - 1), 5);
                EXPECT_EQ(queen_hlooped.neighbors_count((r + 1) * shape[1] - 1), 8);
                EXPECT_EQ(queen_vlooped.neighbors_count((r + 1) * shape[1] - 1), 5);
                EXPECT_EQ(queen_looped.neighbors_count((r + 1) * shape[1] - 1), 8);
            }

            for (std::size_t r = 1; r < shape[0] - 1; ++r)
            {
                for (std::size_t c = 1; c < shape[1] - 1; ++c)
                {
                    // Inner nodes
                    EXPECT_EQ(queen_fixed.neighbors_count(r * shape[1] + c), 8);
                    EXPECT_EQ(queen_hlooped.neighbors_count(r * shape[1] + c), 8);
                    EXPECT_EQ(queen_vlooped.neighbors_count(r * shape[1] + c), 8);
                    EXPECT_EQ(queen_looped.neighbors_count(r * shape[1] + c), 8);
                }
            }
        }
    }
}
