#include "test_raster_grid.hpp"


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace testing
    {

        TEST_F(rook_raster_grid, n_neighbors_max)
        {
            EXPECT_EQ(grid_type::n_neighbors_max(), 4u);
        }

#define EXPECT_INDICES(ROW, COL, INDICES, RC_INDICES)                                              \
    EXPECT_EQ(rook_fixed.neighbors_indices(ROW* shape[1] + COL),                                   \
              (xt::xtensor<std::size_t, 1> INDICES));                                              \
                                                                                                   \
    rook_fixed.neighbors_indices(ROW* shape[1] + COL, neighbors_idx);                              \
    EXPECT_EQ(neighbors_idx, (xt::xtensor<std::size_t, 1> INDICES));                               \
                                                                                                   \
    EXPECT_EQ((rook_fixed.neighbors_indices(ROW, COL)),                                            \
              (grid_type::neighbors_indices_raster_type RC_INDICES));

        TEST_F(rook_raster_grid_fixed, neighbors_indices)
        {
            {  // First row
                // Top-left corner
                EXPECT_INDICES(0, 0, ({ 1, 10 }), ({ /*(r,c)*/ { 0, 1 }, { 1, 0 } }));
                // Inner cols
                for (std::size_t c = 1; c < 9; ++c)
                {
                    EXPECT_INDICES(0,
                                   c,
                                   ({ c - 1, c + 1, shape[1] + c }),
                                   ({
                                       { 0, c - 1 },
                                       /*(r,c)*/ { 0, c + 1 },
                                       { 1, c },
                                   }));
                }
                // Top-right corner
                EXPECT_INDICES(0,
                               9,
                               ({ 8, 19 }),
                               ({ { 0, 8 }, /*(r,c)*/
                                  { 1, 9 } }));
            }
            {  // Inners rows
                for (std::size_t r = 1; r < 4; ++r)
                {
                    // First col
                    EXPECT_INDICES(r,
                                   0,
                                   ({ (r - 1) * shape[1], r * shape[1] + 1, (r + 1) * shape[1] }),
                                   ({ { r - 1, 0 },
                                      /*(r,c)*/ { r, 1 },
                                      { r + 1, 0 } }));
                    // Inners cols
                    for (std::size_t c = 1; c < 9; ++c)
                    {
                        EXPECT_INDICES(r,
                                       c,
                                       ({ (r - 1) * shape[1] + c,
                                          r * shape[1] + c - 1,
                                          r * shape[1] + c + 1,
                                          (r + 1) * shape[1] + c }),
                                       ({
                                           { r - 1, c },
                                           { r, c - 1 },
                                           /*(r,c)*/ { r, c + 1 },
                                           { r + 1, c },
                                       }));
                    }
                    // last col
                    EXPECT_INDICES(
                        r,
                        shape[1] - 1,
                        ({ (r - 1) * shape[1] + 9, r * shape[1] + 8, (r + 1) * shape[1] + 9 }),
                        ({ { r - 1, 9 },
                           { r, 8 }, /*(r,c)*/
                           { r + 1, 9 } }));
                }
            }
            {  // Last row
                // Bottom-left corner
                EXPECT_INDICES(4,
                               0,
                               ({ 30, 41 }),
                               ({ { 3, 0 },
                                  /*(r,c)*/ { 4, 1 } }));
                // Inner cols
                for (std::size_t c = 1; c < 9; ++c)
                {
                    EXPECT_INDICES(4,
                                   c,
                                   ({ 30 + c, 40 + c - 1, 40 + c + 1 }),
                                   ({ { 3, c }, { 4, c - 1 }, /*(r,c)*/ { 4, c + 1 } }));
                }
                // Bottom-right corner
                EXPECT_INDICES(4, 9, ({ 39, 48 }), ({ { 3, 9 }, { 4, 8 } /*(r,c)*/ }));
            }
        }
#undef EXPECT_INDICES

#define EXPECT_INDICES(ROW, COL, INDICES, RC_INDICES)                                              \
    EXPECT_EQ(rook_looped.neighbors_indices(ROW* shape[1] + COL),                                  \
              (xt::xtensor<std::size_t, 1> INDICES));                                              \
                                                                                                   \
    rook_looped.neighbors_indices(ROW* shape[1] + COL, neighbors_idx);                             \
    EXPECT_EQ(neighbors_idx, (xt::xtensor<std::size_t, 1> INDICES));                               \
                                                                                                   \
    EXPECT_EQ((rook_looped.neighbors_indices(ROW, COL)),                                           \
              (grid_type::neighbors_indices_raster_type RC_INDICES));

        TEST_F(rook_raster_grid_looped, neighbors_indices)
        {
            {  // First row
                // Top-left corner
                EXPECT_INDICES(0,
                               0,
                               ({ 40, 9, 1, 10 }),
                               ({ { 4, 0 }, { 0, 9 }, /*(r,c)*/ { 0, 1 }, { 1, 0 } }));
                // Inner cols
                for (std::size_t c = 1; c < 9; ++c)
                {
                    EXPECT_INDICES(0,
                                   c,
                                   ({ 40 + c, c - 1, c + 1, 10 + c }),
                                   ({ { 4, c }, { 0, c - 1 }, /*(r,c)*/ { 0, c + 1 }, { 1, c } }));
                }
                // Top-right corner
                EXPECT_INDICES(0,
                               9,
                               ({ 49, 8, 0, 19 }),
                               ({ { 4, 9 }, { 0, 8 }, /*(r,c)*/ { 0, 0 }, { 1, 9 } }));
            }
            {  // Inners rows
                for (std::size_t r = 1; r < 4; ++r)
                {
                    // First col
                    EXPECT_INDICES(r,
                                   0,
                                   ({ (r - 1) * shape[1],
                                      r * shape[1] + 9,
                                      r * shape[1] + 1,
                                      (r + 1) * shape[1] }),
                                   ({ { r - 1, 0 }, { r, 9 }, /*(r,c)*/ { r, 1 }, { r + 1, 0 } }));
                    // Inners cols
                    for (std::size_t c = 1; c < 9; ++c)
                    {
                        EXPECT_INDICES(r,
                                       c,
                                       ({ (r - 1) * shape[1] + c,
                                          r * shape[1] + c - 1,
                                          r * shape[1] + c + 1,
                                          (r + 1) * shape[1] + c }),
                                       ({
                                           { r - 1, c },
                                           { r, c - 1 },
                                           /*(r,c)*/ { r, c + 1 },
                                           { r + 1, c },
                                       }));
                    }
                    // last col
                    EXPECT_INDICES(r,
                                   shape[1] - 1,
                                   ({ (r - 1) * shape[1] + 9,
                                      r * shape[1] + 8,
                                      r * shape[1],
                                      (r + 1) * shape[1] + 9 }),
                                   ({ { r - 1, 9 }, { r, 8 }, /*(r,c)*/ { r, 0 }, { r + 1, 9 } }));
                }
            }
            {  // Last row
                // Bottom-left corner
                EXPECT_INDICES(4,
                               0,
                               ({ 30, 49, 41, 0 }),
                               ({ { 3, 0 }, { 4, 9 }, /*(r,c)*/ { 4, 1 }, { 0, 0 } }));
                // Inner cols
                for (std::size_t c = 1; c < 9; ++c)
                {
                    EXPECT_INDICES(
                        4,
                        c,
                        ({ 3 * shape[1] + c, 4 * shape[1] + c - 1, 4 * shape[1] + c + 1, c }),
                        ({ { 3, c }, { 4, c - 1 }, /*(r,c)*/ { 4, c + 1 }, { 0, c } }));
                }
                // Bottom-right corner
                EXPECT_INDICES(4,
                               9,
                               ({ 39, 48, 40, 9 }),
                               ({ { 3, 9 }, { 4, 8 }, /*(r,c)*/ { 4, 0 }, { 0, 9 } }));
            }
        }
#undef EXPECT_INDICES

        TEST_F(rook_raster_grid_fixed, neighbors_distances)
        {
            using neighbors_distances_type = grid_type::neighbors_distances_type;

            // Top-left corner
            EXPECT_EQ(rook_fixed.neighbors_distances(0), (neighbors_distances_type{ 1.2, 1.3 }));

            // Bottom-right corner
            EXPECT_EQ(rook_fixed.neighbors_distances(49), (neighbors_distances_type{ 1.3, 1.2 }));

            // Inners rows/cols
            for (std::size_t r = 1; r < 4; ++r)
            {
                for (std::size_t c = 1; c < 9; ++c)
                {
                    EXPECT_EQ(rook_fixed.neighbors_distances(r * shape[1] + c),
                              (neighbors_distances_type{ 1.3, 1.2, 1.2, 1.3 }));
                }
            }
        }

        TEST_F(rook_raster_grid_looped, neighbors_distances)
        {
            using neighbors_distances_type = grid_type::neighbors_distances_type;

            for (std::size_t r = 0; r < 5; ++r)
            {
                for (std::size_t c = 0; c < 10; ++c)
                {
                    EXPECT_EQ(rook_looped.neighbors_distances(r * shape[1] + c),
                              (neighbors_distances_type{ 1.4, 1.8, 1.8, 1.4 }));
                }
            }
        }

        TEST_F(rook_raster_grid_looped, neighbors)
        {  // Do not test extensively, indices and distances are already tested
            using neighbors_type = std::vector<fs::neighbor>;

            {  // First row
                // Top-left corner
                EXPECT_EQ(rook_fixed.neighbors(0),
                          (neighbors_type{
                              /* Node */ { 1, 1.2, fb },
                              { 10, 1.3, fb },
                          }));

                EXPECT_EQ(rook_looped.neighbors(0),
                          (neighbors_type{ { 40, 1.4, lb },
                                           { 9, 1.8, lb },
                                           /* Node */ { 1, 1.8, lb },
                                           { 10, 1.4, lb } }));

                // Inner cols
                for (std::size_t c = 2; c < 8; ++c)
                {
                    EXPECT_EQ(rook_fixed.neighbors(c),
                              (neighbors_type{ { c - 1, 1.2, fb },
                                               /* Node */ { c + 1, 1.2, fb },
                                               { c + 10, 1.3, co } }));

                    EXPECT_EQ(rook_looped.neighbors(c),
                              (neighbors_type{ { c + 40, 1.4, lb },
                                               { c - 1, 1.8, lb },
                                               /* Node */ { c + 1, 1.8, lb },
                                               { c + 10, 1.4, co } }));
                }

                EXPECT_EQ(rook_looped.neighbors(22),
                          (neighbors_type{ { 12, 1.4, co },
                                           { 21, 1.8, co },
                                           /* Node */ { 23, 1.8, co },
                                           { 32, 1.4, co } }));
            }
        }

        TEST_F(rook_raster_grid, neighbors_count)
        {
            // Top-left corner nodes
            EXPECT_EQ(rook_fixed.neighbors_count(static_cast<std::size_t>(0)), 2);
            EXPECT_EQ(rook_hlooped.neighbors_count(static_cast<std::size_t>(0)), 3);
            EXPECT_EQ(rook_vlooped.neighbors_count(static_cast<std::size_t>(0)), 3);
            EXPECT_EQ(rook_looped.neighbors_count(static_cast<std::size_t>(0)), 4);

            // Top-right corner nodes
            EXPECT_EQ(rook_fixed.neighbors_count(shape[1] - 1), 2);
            EXPECT_EQ(rook_hlooped.neighbors_count(shape[1] - 1), 3);
            EXPECT_EQ(rook_vlooped.neighbors_count(shape[1] - 1), 3);
            EXPECT_EQ(rook_looped.neighbors_count(shape[1] - 1), 4);

            // Bottom-left corner nodes
            EXPECT_EQ(rook_fixed.neighbors_count((shape[0] - 1) * shape[1]), 2);
            EXPECT_EQ(rook_hlooped.neighbors_count((shape[0] - 1) * shape[1]), 3);
            EXPECT_EQ(rook_vlooped.neighbors_count((shape[0] - 1) * shape[1]), 3);
            EXPECT_EQ(rook_looped.neighbors_count((shape[0] - 1) * shape[1]), 4);

            // Bottom-right corner nodes
            EXPECT_EQ(rook_fixed.neighbors_count(shape[0] * shape[1] - 1), 2);
            EXPECT_EQ(rook_hlooped.neighbors_count(shape[0] * shape[1] - 1), 3);
            EXPECT_EQ(rook_vlooped.neighbors_count(shape[0] * shape[1] - 1), 3);
            EXPECT_EQ(rook_looped.neighbors_count(shape[0] * shape[1] - 1), 4);

            for (std::size_t c = 1; c < shape[1] - 1; ++c)
            {
                // Top edge nodes (without corners)
                EXPECT_EQ(rook_fixed.neighbors_count(c), 3);
                EXPECT_EQ(rook_hlooped.neighbors_count(c), 3);
                EXPECT_EQ(rook_vlooped.neighbors_count(c), 4);
                EXPECT_EQ(rook_looped.neighbors_count(c), 4);

                // Bottom edge nodes (without corners)
                EXPECT_EQ(rook_fixed.neighbors_count((shape[0] - 1) * shape[1] + c), 3);
                EXPECT_EQ(rook_hlooped.neighbors_count((shape[0] - 1) * shape[1] + c), 3);
                EXPECT_EQ(rook_vlooped.neighbors_count((shape[0] - 1) * shape[1] + c), 4);
                EXPECT_EQ(rook_looped.neighbors_count((shape[0] - 1) * shape[1] + c), 4);
            }

            for (std::size_t r = 1; r < shape[0] - 1; ++r)
            {
                // Left edge nodes (without corners)
                EXPECT_EQ(rook_fixed.neighbors_count(r * shape[1]), 3);
                EXPECT_EQ(rook_hlooped.neighbors_count(r * shape[1]), 4);
                EXPECT_EQ(rook_vlooped.neighbors_count(r * shape[1]), 3);
                EXPECT_EQ(rook_looped.neighbors_count(r * shape[1]), 4);

                // Right edge nodes (without corners)
                EXPECT_EQ(rook_fixed.neighbors_count((r + 1) * shape[1] - 1), 3);
                EXPECT_EQ(rook_hlooped.neighbors_count((r + 1) * shape[1] - 1), 4);
                EXPECT_EQ(rook_vlooped.neighbors_count((r + 1) * shape[1] - 1), 3);
                EXPECT_EQ(rook_looped.neighbors_count((r + 1) * shape[1] - 1), 4);
            }

            for (std::size_t r = 1; r < shape[0] - 1; ++r)
            {
                for (std::size_t c = 1; c < shape[1] - 1; ++c)
                {
                    // Inner nodes
                    EXPECT_EQ(rook_fixed.neighbors_count(r * shape[1] + c), 4);
                    EXPECT_EQ(rook_hlooped.neighbors_count(r * shape[1] + c), 4);
                    EXPECT_EQ(rook_vlooped.neighbors_count(r * shape[1] + c), 4);
                    EXPECT_EQ(rook_looped.neighbors_count(r * shape[1] + c), 4);
                }
            }
        }
    }
}
