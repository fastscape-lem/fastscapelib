#include "xtensor/xtensor.hpp"

#include "fastscapelib/algo/pflood.hpp"

#include "test_pflood.hpp"


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace testing
    {
        TEST_F(sinks_raster_grid, fill_sinks_sloped)
        {
            using elev_type = xt::xtensor<double, 2>;

            {
                SCOPED_TRACE("full closed boundaries");

                elev_type res = elevation;
                fs::detail::fill_sinks_sloped(graph_full_closed.impl(), res);

                // center node "filled"
                EXPECT_GT(res(1, 1), res(0, 1));
                // all border nodes should have elevation unchanged
                EXPECT_EQ(res(0, 1), elevation(0, 1));
            }

            {
                SCOPED_TRACE("left closed boundary");

                elev_type res = elevation;
                fs::detail::fill_sinks_sloped(graph_left_closed.impl(), res);

                // center node "filled"
                EXPECT_GT(res(1, 1), res(0, 0));
                // mid-top node "filled"
                EXPECT_GT(res(0, 1), res(0, 0));
            }

            {
                SCOPED_TRACE("vertical looped boundaries");

                elev_type res = elevation;
                fs::detail::fill_sinks_sloped(graph_vert_looped.impl(), res);

                // center node "filled"
                EXPECT_GT(res(1, 1), res(0, 0));
                // mid-top node "filled"
                EXPECT_GT(res(0, 1), res(0, 0));
            }
        }

        TEST_F(sinks_profile_grid, fill_sinks_sloped)
        {
            using elev_type = xt::xtensor<double, 1>;

            {
                SCOPED_TRACE("closed boundaries");

                elev_type res = elevation;
                fs::detail::fill_sinks_sloped(graph_closed.impl(), res);

                EXPECT_EQ(res(0), elevation(0));
                EXPECT_EQ(res(3), elevation(3));
                EXPECT_GT(res(2), res(3));
                EXPECT_GT(res(1), res(2));
            }

            {
                SCOPED_TRACE("half-open boundaries");

                elev_type res = elevation;
                fs::detail::fill_sinks_sloped(graph_half_open.impl(), res);

                EXPECT_EQ(res(0), elevation(0));
                EXPECT_GT(res(1), res(0));
                EXPECT_GT(res(2), res(1));
                EXPECT_GT(res(3), res(2));
            }

            {
                SCOPED_TRACE("closed boundaries and mask");

                xt::xtensor<bool, 1> mask{ false, true, false, false };
                graph_closed.set_mask(mask);
                elev_type res = elevation;
                fs::detail::fill_sinks_sloped(graph_closed.impl(), res);

                EXPECT_EQ(res(0), elevation(0));
                // masked node not filled
                EXPECT_EQ(res(1), elevation(1));
                EXPECT_EQ(res(3), elevation(3));
                EXPECT_GT(res(2), res(3));
                // masked node not filled
                EXPECT_GT(res(2), res(1));
            }
        }
    }
}
