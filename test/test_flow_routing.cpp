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

    xt::xtensor<index_t, 1> receivers = xt::ones<index_t>({16}) * -1;

    xt::xtensor<index_t, 1> expected_receivers
        { 0,  1,  2,  3,
          4,  9,  7,  7,
          8,  9,  9, 11,
         12, 13, 14, 15};

    xt::xtensor<double, 1> dist2receivers = xt::ones<double>({16}) * -1.;

    xt::xtensor<double, 1> expected_dist2receivers
        {0.,  0.,  0.,  0.,
         0.,  1.,  1.,  0.,
         0.,  0.,  1.,  0.,
         0.,  0.,  0.,  0.};

    fs::compute_receivers_d8(receivers, dist2receivers,
                             elevation, active_nodes,
                             1., 1.);

    //TODO: test with different dx, dy
    //TODO: test case with diagonal receivers (check dist2receivers)
    //TODO: -> ideally, test case should include all 8 directions
    //      (maybe create fixtures with transpose)

    ASSERT_TRUE(xt::all(xt::equal(receivers, expected_receivers)));
    ASSERT_TRUE(xt::allclose(dist2receivers, expected_dist2receivers));
}


TEST(flow_routing, compute_donors) {
    // Example in Braun and Willet, 2013 as a test case.
    xt::xtensor<index_t, 1> receivers {1, 4, 1, 6, 4, 4, 5, 4, 6, 7};
    xt::xtensor<index_t, 1> ndonors = xt::ones<index_t>({10}) * -1;
    xt::xtensor<index_t, 2> donors = xt::ones<index_t>({10, 8}) * -1;

    xt::xtensor<index_t, 1> expected_ndonors {0, 2, 0, 0, 3, 1, 2, 1, 0, 0};
    xt::xtensor<index_t, 2> expected_donors
        {{-1, -1, -1, -1, -1, -1, -1, -1},
         { 0,  2, -1, -1, -1, -1, -1, -1},
         {-1, -1, -1, -1, -1, -1, -1, -1},
         {-1, -1, -1, -1, -1, -1, -1, -1},
         { 1,  5,  7, -1, -1, -1, -1, -1},
         { 6, -1, -1, -1, -1, -1, -1, -1},
         { 3,  8, -1, -1, -1, -1, -1, -1},
         { 9, -1, -1, -1, -1, -1, -1, -1},
         {-1, -1, -1, -1, -1, -1, -1, -1},
         {-1, -1, -1, -1, -1, -1, -1, -1}};

    fs::compute_donors(ndonors, donors, receivers);

    ASSERT_TRUE(xt::all(xt::equal(ndonors, expected_ndonors)));
    ASSERT_TRUE(xt::all(xt::equal(donors, expected_donors)));
}


TEST(flow_routing, compute_stack) {
    // Example in Braun and Willet, 2013 as a test case.
    xt::xtensor<index_t, 1> receivers {1, 4, 1, 6, 4, 4, 5, 4, 6, 7};
    xt::xtensor<index_t, 1> ndonors   {0, 2, 0, 0, 3, 1, 2, 1, 0, 0};
    xt::xtensor<index_t, 2> donors
        {{-1, -1, -1, -1, -1, -1, -1, -1},
         { 0,  2, -1, -1, -1, -1, -1, -1},
         {-1, -1, -1, -1, -1, -1, -1, -1},
         {-1, -1, -1, -1, -1, -1, -1, -1},
         { 1,  5,  7, -1, -1, -1, -1, -1},
         { 6, -1, -1, -1, -1, -1, -1, -1},
         { 3,  8, -1, -1, -1, -1, -1, -1},
         { 9, -1, -1, -1, -1, -1, -1, -1},
         {-1, -1, -1, -1, -1, -1, -1, -1},
         {-1, -1, -1, -1, -1, -1, -1, -1}};

    xt::xtensor<index_t, 1> stack = xt::ones<index_t>({10}) * -1;
    xt::xtensor<index_t, 1> expected_stack {4, 1, 0, 2, 5, 6, 3, 8, 7, 9};

    fs::compute_stack(stack, ndonors, donors, receivers);

    ASSERT_TRUE(xt::all(xt::equal(stack, expected_stack)));
}
