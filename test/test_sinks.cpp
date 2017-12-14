#include "gtest/gtest.h"
#include "xtensor/xtensor.hpp"

#include "fastscapelib/sinks.hpp"


TEST(sinks, fill_sinks_flat)
{
    xt::xtensor<double, 2> arr
       {{1.0, 2.0, 3.0},
        {2.0, 0.1, 7.0},
        {2.0, 5.0, 7.0}};

    fs::fill_sinks_flat(arr);

    EXPECT_EQ(arr(1, 1), arr(0, 0));
}


TEST(sinks, fill_sinks_sloped)
{
    xt::xtensor<double, 2> arr
       {{1.0, 2.0, 3.0},
        {2.0, 0.1, 7.0},
        {2.0, 5.0, 7.0}};

    fs::fill_sinks_sloped(arr);

    EXPECT_GT(arr(1, 1), arr(0, 0));
}
