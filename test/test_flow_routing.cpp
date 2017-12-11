#include "gtest/gtest.h"
#include "xtensor/xtensor.hpp"

#include "fastscape/utils.hpp"
#include "fastscape/flow_routing.hpp"


TEST(flow_routing, compute_receivers_d8) {
    xt::xtensor<double, 2> elevation
        {{0.82,  0.16,  0.14,  0.20},
         {0.71,  0.97,  0.41,  0.09},
         {0.49,  0.01,  0.19,  0.38},
         {0.29,  0.82,  0.09,  0.88}};

    xt::xtensor<bool, 2> active_nodes
       {{false,  false,  false,  false},
        {false,  true,   true,   false},
        {false,  true,   true,   false},
        {false,  false,  false,  false}};

    xt::xtensor<index_t, 1> receivers = xt::ones<index_t>({16});

    xt::xtensor<index_t, 1> expected_receivers
        { 0,  1,  2,  3,
          4,  9,  7,  7,
          8,  9,  9, 11,
         12, 13, 14, 15};

    xt::xtensor<double, 1> dist2receivers = xt::zeros<double>({16});

    xt::xtensor<double, 1> expected_dist2receivers
        {0.,  0.,  0.,  0.,
         0.,  1.,  1.,  0.,
         0.,  0.,  1.,  0.,
         0.,  0.,  0.,  0.};

    fs::compute_receivers_d8(receivers, dist2receivers,
                             elevation, active_nodes,
                             1., 1.);

    //TODO: test with different dx, dy

    ASSERT_TRUE(xt::all(xt::equal(receivers, expected_receivers)));
    ASSERT_TRUE(xt::allclose(dist2receivers, expected_dist2receivers));
}
