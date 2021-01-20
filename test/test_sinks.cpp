#include "xtensor/xtensor.hpp"

#include "fastscapelib/sinks.hpp"

#include "test_sinks.hpp"


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace testing
    {
        TEST_F(sinks_raster_grid, fill_sinks_flat)
        {
            using elev_type = xt::xtensor<double, 2>;

            {
                SCOPED_TRACE("full closed boundaries");

                elev_type res = elevation;
                fs::fill_sinks_flat(raster_grid_full_closed, res);

                // center node "filled"
                EXPECT_EQ(res(1, 1), res(0, 1));
                // all border nodes should have elevation unchanged
                EXPECT_EQ(res(0, 1), elevation(0, 1));
            }

            {
                SCOPED_TRACE("left closed boundary");

                elev_type res = elevation;
                fs::fill_sinks_flat(raster_grid_left_closed, res);

                // center node "filled"
                EXPECT_EQ(res(1, 1), res(0, 0));
                // mid-top node "filled"
                EXPECT_EQ(res(0, 1), res(0, 0));
            }

            {
                SCOPED_TRACE("vertical looped boundaries");

                elev_type res = elevation;
                fs::fill_sinks_flat(raster_grid_vert_looped, res);

                // center node "filled"
                EXPECT_EQ(res(1, 1), res(0, 0));
                // mid-top node "filled"
                EXPECT_EQ(res(0, 1), res(0, 0));
            }
        }

        TEST_F(sinks_raster_grid, fill_sinks_sloped)
        {
            using elev_type = xt::xtensor<double, 2>;

            {
                SCOPED_TRACE("full closed boundaries");

                elev_type res = elevation;
                fs::fill_sinks_sloped(raster_grid_full_closed, res);

                // center node "filled"
                EXPECT_GT(res(1, 1), res(0, 1));
                // all border nodes should have elevation unchanged
                EXPECT_EQ(res(0, 1), elevation(0, 1));
            }

            {
                SCOPED_TRACE("left closed boundary");

                elev_type res = elevation;
                fs::fill_sinks_sloped(raster_grid_left_closed, res);

                // center node "filled"
                EXPECT_GT(res(1, 1), res(0, 0));
                // mid-top node "filled"
                EXPECT_GT(res(0, 1), res(0, 0));
            }

            {
                SCOPED_TRACE("vertical looped boundaries");

                elev_type res = elevation;
                fs::fill_sinks_sloped(raster_grid_vert_looped, res);

                // center node "filled"
                EXPECT_GT(res(1, 1), res(0, 0));
                // mid-top node "filled"
                EXPECT_GT(res(0, 1), res(0, 0));
            }
        }

        TEST_F(sinks_profile_grid, fill_sinks_flat)
        {
            using elev_type = xt::xtensor<double, 1>;

            {
                SCOPED_TRACE("closed boundaries");

                elev_type res = elevation;
                fs::fill_sinks_flat(profile_grid_closed, res);

                elev_type expected {3.0, 2.0, 2.0, 2.0};
                EXPECT_TRUE(xt::all(xt::equal(res, expected)));
            }

            {
                SCOPED_TRACE("half-open boundaries");

                elev_type res = elevation;
                fs::fill_sinks_flat(profile_grid_half_open, res);

                elev_type expected {3.0, 3.0, 3.0, 3.0};
                EXPECT_TRUE(xt::all(xt::equal(res, expected)));
            }
        }

        TEST_F(sinks_profile_grid, fill_sinks_sloped)
        {
            using elev_type = xt::xtensor<double, 1>;

            {
                SCOPED_TRACE("closed boundaries");

                elev_type res = elevation;
                fs::fill_sinks_sloped(profile_grid_closed, res);

                EXPECT_EQ(res(0), elevation(0));
                EXPECT_EQ(res(3), elevation(3));
                EXPECT_GT(res(2), res(3));
                EXPECT_GT(res(1), res(2));
            }

            {
                SCOPED_TRACE("half-open boundaries");

                elev_type res = elevation;
                fs::fill_sinks_sloped(profile_grid_half_open, res);

                EXPECT_EQ(res(0), elevation(0));
                EXPECT_GT(res(1), res(0));
                EXPECT_GT(res(2), res(1));
                EXPECT_GT(res(3), res(2));
            }
        }
    }
}
