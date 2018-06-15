#include "gtest/gtest.h"
#include "xtensor/xtensor.hpp"

#include "fastscapelib/bedrock_channel.hpp"


namespace fs = fastscapelib;


TEST(bedrock_channel, erode_stream_power)
{
    // Example in Braun and Willet, 2013 as a test case
    // - fixed drainage area = 2 for all cells
    // - fixed dist2receivers = 1 for all nodes (except outlet)
    xt::xtensor<index_t, 1> receivers {1, 4, 1, 6, 4, 4, 5, 4, 6, 7};
    xt::xtensor<double, 1> dist2receivers = xt::ones<double>({10}) * 2.;
    xt::xtensor<index_t, 1> stack {4, 1, 0, 2, 5, 6, 3, 8, 7, 9};
    xt::xtensor<double, 1> drainage_area = xt::ones<double>({10}) * 2.;
    xt::xtensor<double, 1> elevation = xt::zeros<double>({10});
    xt::xtensor<double, 1> erosion = xt::zeros<double>({10});

    fs::erode_stream_power(erosion, elevation, stack,
                           receivers, dist2receivers, drainage_area,
                           1e-7, 0.4, 1., 1000., 1e-3);

    //EXPECT_EQ(arr(1, 1), arr(0, 0));
}
