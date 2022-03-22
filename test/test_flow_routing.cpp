#include "fastscapelib/utils/utils.hpp"
#include "fastscapelib/grid/raster_grid.hpp"
#include "fastscapelib/flow_routing.hpp"

#include "xtensor/xtensor.hpp"
#include "xtensor/xio.hpp"

#include "gtest/gtest.h"

#include <cmath>


namespace fs = fastscapelib;


TEST(flow_routing, get_d8_distances)
{
    auto d8_dist = fs::detail::get_d8_distances(2., 1.);

    double ddiag = std::sqrt(4 + 1);

    std::array<double, 9> expected{ 0., 1., ddiag, 2., ddiag, 1., ddiag, 2., ddiag };

    EXPECT_TRUE(std::equal(d8_dist.begin(), d8_dist.end(), expected.begin(), expected.end()));

    // TODO: consider using google mock which gives more details about failure
    // EXPECT_THAT(d8_dist, ::testing::ContainerEq(expected));
}


TEST(flow_routing, compute_receivers_d8)
{
    xt::xtensor<double, 2> elevation{ { 0.82, 0.16, 0.14, 0.20 },
                                      { 0.71, 0.97, 0.41, 0.09 },
                                      { 0.49, 0.01, 0.19, 0.38 },
                                      { 0.29, 0.82, 0.09, 0.88 } };

    xt::xtensor<bool, 2> active_nodes{ { false, false, false, false },
                                       { false, true, true, false },
                                       { false, true, true, false },
                                       { false, false, false, false } };

    xt::xtensor<index_t, 1> receivers = xt::ones<index_t>({ 16 }) * -1;

    xt::xtensor<index_t, 1> expected_receivers{
        0, 1, 2, 3, 4, 9, 7, 7, 8, 9, 9, 11, 12, 13, 14, 15
    };

    xt::xtensor<double, 1> dist2receivers = xt::ones<double>({ 16 }) * -1.;

    xt::xtensor<double, 1> expected_dist2receivers{ 0., 0., 0., 0., 0., 1., 1., 0.,
                                                    0., 0., 1., 0., 0., 0., 0., 0. };

    fs::compute_receivers_d8(receivers, dist2receivers, elevation, active_nodes, 1., 1.);

    // TODO: test with different dx, dy
    // TODO: test case with diagonal receivers (check dist2receivers)
    // TODO: -> ideally, test case should include all 8 directions
    //       (maybe create fixtures with transpose)

    EXPECT_TRUE(xt::all(xt::equal(receivers, expected_receivers)));
    EXPECT_TRUE(xt::allclose(dist2receivers, expected_dist2receivers));
}


TEST(flow_routing, compute_donors)
{
    // Example in Braun and Willet, 2013 as a test case.
    xt::xtensor<index_t, 1> receivers{ 1, 4, 1, 6, 4, 4, 5, 4, 6, 7 };
    xt::xtensor<index_t, 1> ndonors = xt::ones<index_t>({ 10 }) * -1;
    xt::xtensor<index_t, 2> donors = xt::ones<index_t>({ 10, 8 }) * -1;

    xt::xtensor<index_t, 1> expected_ndonors{ 0, 2, 0, 0, 3, 1, 2, 1, 0, 0 };
    xt::xtensor<index_t, 2> expected_donors{
        { -1, -1, -1, -1, -1, -1, -1, -1 }, { 0, 2, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1, -1, -1, -1 },
        { 1, 5, 7, -1, -1, -1, -1, -1 },    { 6, -1, -1, -1, -1, -1, -1, -1 },
        { 3, 8, -1, -1, -1, -1, -1, -1 },   { 9, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1, -1, -1, -1 }
    };

    fs::compute_donors(ndonors, donors, receivers);

    EXPECT_TRUE(xt::all(xt::equal(ndonors, expected_ndonors)));
    EXPECT_TRUE(xt::all(xt::equal(donors, expected_donors)));
}


TEST(flow_routing, compute_stack)
{
    // Example in Braun and Willet, 2013 as a test case.
    xt::xtensor<index_t, 1> receivers{ 1, 4, 1, 6, 4, 4, 5, 4, 6, 7 };
    xt::xtensor<index_t, 1> ndonors{ 0, 2, 0, 0, 3, 1, 2, 1, 0, 0 };
    xt::xtensor<index_t, 2> donors{
        { -1, -1, -1, -1, -1, -1, -1, -1 }, { 0, 2, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1, -1, -1, -1 },
        { 1, 5, 7, -1, -1, -1, -1, -1 },    { 6, -1, -1, -1, -1, -1, -1, -1 },
        { 3, 8, -1, -1, -1, -1, -1, -1 },   { 9, -1, -1, -1, -1, -1, -1, -1 },
        { -1, -1, -1, -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1, -1, -1, -1 }
    };

    xt::xtensor<index_t, 1> stack = xt::ones<index_t>({ 10 }) * -1;
    xt::xtensor<index_t, 1> expected_stack{ 4, 1, 0, 2, 5, 6, 3, 8, 7, 9 };

    fs::compute_stack(stack, ndonors, donors, receivers);

    EXPECT_TRUE(xt::all(xt::equal(stack, expected_stack)));
}


TEST(flow_routing, compute_basins)
{
    // Example in Braun and Willet, 2013 as a test case.
    // TODO: test with multiple basins in node network.
    xt::xtensor<index_t, 1> receivers{ 1, 4, 1, 6, 4, 4, 5, 4, 6, 7 };
    xt::xtensor<index_t, 1> stack{ 4, 1, 0, 2, 5, 6, 3, 8, 7, 9 };
    xt::xtensor<index_t, 1> basins = xt::ones<index_t>({ 10 }) * -1;
    xt::xtensor<index_t, 1> outlets = xt::ones<index_t>({ 10 }) * -1;

    xt::xtensor<index_t, 1> expected_basins{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    xt::xtensor<index_t, 1> expected_outlets{ 4, -1, -1, -1, -1, -1, -1, -1, -1, -1 };

    index_t nbasins = fs::compute_basins(basins, outlets, stack, receivers);

    EXPECT_EQ(nbasins, 1);
    EXPECT_TRUE(xt::all(xt::equal(basins, expected_basins)));
    EXPECT_TRUE(xt::all(xt::equal(outlets, expected_outlets)));
}


TEST(flow_routing, find_pits)
{
    // simple 4x4 test case with fixed boundaries and pit located at (2, 1)
    // 12 boundary nodes (i.e., 1-node open basin) + 1 pit = 13 basins
    index_t nbasins = 13;
    xt::xtensor<index_t, 1> outlets{ 0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15, 9, -1, -1, -1 };

    xt::xtensor<bool, 2> active_nodes{ { false, false, false, false },
                                       { false, true, true, false },
                                       { false, true, true, false },
                                       { false, false, false, false } };

    xt::xtensor<index_t, 1> pits = xt::ones<index_t>({ 16 }) * -1;

    xt::xtensor<index_t, 1> expected_pits{ 9,  -1, -1, -1, -1, -1, -1, -1,
                                           -1, -1, -1, -1, -1, -1, -1, -1 };

    index_t npits = fs::find_pits(pits, outlets, active_nodes, nbasins);

    EXPECT_EQ(npits, 1);
    EXPECT_TRUE(xt::all(xt::equal(pits, expected_pits)));
}


TEST(flow_routing, compute_drainage_area)
{
    // Example in Braun and Willet, 2013 as a test case
    // (setting constant area = 2 for all cells).
    xt::xtensor<index_t, 1> receivers{ 1, 4, 1, 6, 4, 4, 5, 4, 6, 7 };
    xt::xtensor<index_t, 1> stack{ 4, 1, 0, 2, 5, 6, 3, 8, 7, 9 };

    xt::xtensor<double, 1> expected_area{ 2., 6., 2., 2., 20., 8., 6., 4., 2., 2. };

    xt::xtensor<double, 2> expected_area_2d{ { 2., 6., 2., 2., 20. }, { 8., 6., 4., 2., 2. } };

    {
        SCOPED_TRACE("generic mesh");
        xt::xtensor<double, 1> drainage_area = xt::ones<double>({ 10 }) * -1;
        xt::xtensor<double, 1> cell_area = xt::ones<double>({ 10 }) * 2;
        fs::compute_drainage_area(drainage_area, cell_area, stack, receivers);

        EXPECT_TRUE(xt::all(xt::equal(drainage_area, expected_area)));
    }

    {
        SCOPED_TRACE("2d grid flattened");
        xt::xtensor<double, 1> drainage_area = xt::ones<double>({ 10 }) * -1;
        fs::compute_drainage_area(drainage_area, stack, receivers, 1., 2.);

        EXPECT_TRUE(xt::all(xt::equal(drainage_area, expected_area)));
    }

    {
        SCOPED_TRACE("2d grid");
        xt::xtensor<double, 2> drainage_area = xt::ones<double>({ 2, 5 }) * -1;
        fs::compute_drainage_area(drainage_area, stack, receivers, 1., 2.);

        EXPECT_TRUE(xt::all(xt::equal(drainage_area, expected_area_2d)));
    }

    {
        SCOPED_TRACE("2d grid -- 2d cell_area");
        xt::xtensor<double, 2> drainage_area = xt::ones<double>({ 2, 5 }) * -1;
        xt::xtensor<double, 2> cell_area = xt::ones<double>({ 2, 5 }) * 2;
        fs::compute_drainage_area(drainage_area, cell_area, stack, receivers);

        EXPECT_TRUE(xt::all(xt::equal(drainage_area, expected_area_2d)));
    }
}
