#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/grid/healpix_grid.hpp"
#include "fastscapelib/utils/consts.hpp"

#include "xtensor/xtensor.hpp"
#include "xtensor/xmath.hpp"

#include "gtest/gtest.h"


namespace fs = fastscapelib;


namespace fastscapelib
{
    namespace testing
    {

        class healpix : public ::testing::Test
        {
        protected:
            int nside = 32;
            std::size_t size = 12288;

            xt::xtensor<fs::node_status, 1> nodes_status = xt::zeros<fs::node_status>({ size });
        };

        TEST_F(healpix, static_expr)
        {
            EXPECT_EQ(fs::healpix_grid<>::is_structured(), false);
            EXPECT_EQ(fs::healpix_grid<>::is_uniform(), false);
            EXPECT_EQ(fs::healpix_grid<>::n_neighbors_max(), 8u);
            EXPECT_EQ(fs::healpix_grid<>::container_ndims(), 1);
        }

        TEST_F(healpix, ctor)
        {
            auto grid = fs::healpix_grid<>(nside, nodes_status);

            EXPECT_EQ(grid.nside(), nside);
            EXPECT_EQ(grid.size(), size);
            EXPECT_EQ(grid.radius(), fs::numeric_constants<>::EARTH_RADIUS);
        }

        TEST_F(healpix, nodes_areas)
        {
            auto grid = fs::healpix_grid<>(nside, nodes_status);

            double expected = 41509152987.45021;
            EXPECT_NEAR(grid.nodes_areas(0), expected, 1e-5);

            xt::xtensor<double, 1> expected_arr({ grid.size() }, expected);
            EXPECT_TRUE(xt::allclose(grid.nodes_areas(), expected_arr));
        }
    }
}
