#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/grid/healpix_grid.hpp"
#include "fastscapelib/utils/consts.hpp"

#include "xtensor/containers/xtensor.hpp"
#include "xtensor/core/xmath.hpp"

#include "gtest/gtest.h"


namespace fs = fastscapelib;


namespace fastscapelib
{
    namespace testing
    {

        class healpix_grid : public ::testing::Test
        {
        protected:
            int nside = 32;
            std::size_t size = 12288;

            xt::xtensor<fs::node_status, 1> nodes_status = xt::zeros<fs::node_status>({ size });
        };

        TEST_F(healpix_grid, static_expr)
        {
            EXPECT_EQ(fs::healpix_grid<>::is_structured(), false);
            EXPECT_EQ(fs::healpix_grid<>::is_uniform(), false);
            EXPECT_EQ(fs::healpix_grid<>::n_neighbors_max(), 8u);
            EXPECT_EQ(fs::healpix_grid<>::container_ndims(), 1);
        }

        TEST_F(healpix_grid, ctor)
        {
            auto grid = fs::healpix_grid<>(nside, nodes_status);

            EXPECT_EQ(grid.nside(), nside);
            EXPECT_EQ(grid.size(), size);
            EXPECT_EQ(grid.radius(), fs::numeric_constants<>::EARTH_RADIUS_METERS);
        }

        TEST_F(healpix_grid, nodes_status)
        {
            auto grid = fs::healpix_grid<>(nside, nodes_status);

            EXPECT_EQ(grid.nodes_status(), nodes_status);

            {
                SCOPED_TRACE("setting ghost node and updating node status should update neighbors");

                auto nodes_status2 = nodes_status;
                nodes_status2(1) = fs::node_status::ghost;
                auto grid2 = fs::healpix_grid<>(nside, nodes_status2);
                EXPECT_EQ(grid2.neighbors_count(0), 6);

                auto nodes_status3 = nodes_status;
                nodes_status3(0) = fs::node_status::ghost;
                auto grid3 = fs::healpix_grid<>(nside, nodes_status3);
                EXPECT_EQ(grid3.neighbors_count(0), 0);
            }

            {
                SCOPED_TRACE("looped boundary conditions not supported");

                auto nodes_status2 = nodes_status;
                nodes_status2(0) = fs::node_status::looped;
                EXPECT_THROW(fs::healpix_grid<>(nside, nodes_status2), std::invalid_argument);
            }
        }

        TEST_F(healpix_grid, nodes_areas)
        {
            auto grid = fs::healpix_grid<>(nside, nodes_status);

            double expected = 41509152987.45021;
            EXPECT_NEAR(grid.nodes_areas(0), expected, 1e-5);

            xt::xtensor<double, 1> expected_arr({ grid.size() }, expected);
            EXPECT_TRUE(xt::allclose(grid.nodes_areas(), expected_arr));
        }

        TEST_F(healpix_grid, nodes_lonlat)
        {
            auto grid = fs::healpix_grid<>(nside, nodes_status);

            auto actual = grid.nodes_lonlat(0);
            auto actual_arr = grid.nodes_lonlat();
            auto actual_arr0
                = std::make_pair<double&, double&>(actual_arr.first(0), actual_arr.second(0));

            auto expected = std::make_pair<double, double>(0.7853981633974483, 1.5452801164374776);

            EXPECT_EQ(actual, actual_arr0);
            EXPECT_EQ(actual, expected);
        }

        TEST_F(healpix_grid, nodes_xyz)
        {
            auto grid = fs::healpix_grid<>(nside, nodes_status);

            auto actual = grid.nodes_xyz(0);
            auto actual_arr = grid.nodes_xyz();
            auto actual_arr0 = std::make_tuple<double&, double&, double&>(
                std::get<0>(actual_arr)(0), std::get<1>(actual_arr)(0), std::get<2>(actual_arr)(0));

            auto expected = std::make_tuple<double, double, double>(
                114937.4753788243, 114937.4753788243, 6368926.106770833);

            EXPECT_EQ(actual, actual_arr0);
            EXPECT_EQ(actual, expected);
        }

        TEST_F(healpix_grid, neighbors_count)
        {
            auto grid = fs::healpix_grid<>(nside, nodes_status);

            EXPECT_EQ(grid.neighbors_count(0), 7);
            EXPECT_EQ(grid.neighbors_count(1984), 6);
        }

        TEST_F(healpix_grid, neighbors_indices)
        {
            auto grid = fs::healpix_grid<>(nside, nodes_status);

            EXPECT_EQ(grid.neighbors_indices(0),
                      (xt::xtensor<std::size_t, 1>{ 11, 3, 2, 1, 6, 5, 13 }));
        }

        TEST_F(healpix_grid, nodes_distances)
        {
            auto grid = fs::healpix_grid<>(nside, nodes_status);

            xt::xtensor<double, 1> expected({ 0.04752047,
                                              0.03608146,
                                              0.05102688,
                                              0.03608146,
                                              0.04752047,
                                              0.02914453,
                                              0.0510435 });
            EXPECT_TRUE(xt::allclose(grid.neighbors_distances(0), expected));
        }
    }
}
