#include "fastscapelib/grid/healpix_grid.hpp"

#include "fastscapelib/utils/consts.hpp"
#include "xtensor/xtensor.hpp"
#include "gtest/gtest.h"
#include <xtensor/xmath.hpp>


namespace fs = fastscapelib;


namespace fastscapelib
{
    namespace testing
    {

        TEST(healpix_grid, ctor)
        {
            auto grid = fs::healpix_grid<>(32);

            EXPECT_EQ(grid.nside(), 32);
            EXPECT_EQ(grid.size(), 12288);
            EXPECT_EQ(grid.radius(), fs::numeric_constants<>::EARTH_RADIUS);
        }

        TEST(healpix_grid, nodes_areas)
        {
            auto grid = fs::healpix_grid<>(32);

            double expected = 41509152987.45021;
            EXPECT_NEAR(grid.nodes_areas(0), expected, 1e-5);

            xt::xtensor<double, 1> expected_arr({ grid.size() }, expected);
            EXPECT_TRUE(xt::allclose(grid.nodes_areas(), expected_arr));
        }
    }
}
