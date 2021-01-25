#ifndef FASTSCAPELIB_TEST_SINKS_H
#define FASTSCAPELIB_TEST_SINKS_H

#include "xtensor/xtensor.hpp"

#include "fastscapelib/profile_grid.hpp"
#include "fastscapelib/raster_grid.hpp"

#include "gtest/gtest.h"

#include <array>


namespace fs = fastscapelib;


namespace fastscapelib
{
    namespace testing
    {

        class sinks_raster_grid: public ::testing::Test
        {
        protected:

            using grid_type = fs::raster_grid;
            using shape_type = typename grid_type::shape_type;
            using elev_type = xt::xtensor<double, 2>;

            fs::node_status ns_core = fs::node_status::core;
            fs::node_status ns_fixed = fs::node_status::fixed_value_boundary;
            fs::node_status ns_looped = fs::node_status::looped_boundary;

            std::array<fs::node_status, 4> lclosed = {{ ns_fixed, ns_core, ns_core, ns_core }};
            std::array<fs::node_status, 4> vlooped = {{ ns_fixed, ns_core, ns_looped, ns_looped }};

            fs::raster_boundary_status full_closed_status {ns_fixed};
            fs::raster_boundary_status left_closed_status {lclosed};
            fs::raster_boundary_status vert_looped_status {vlooped};

            shape_type shape {{3, 3}};

            grid_type raster_grid_full_closed = grid_type(shape, {1.0, 1.0}, full_closed_status);
            grid_type raster_grid_left_closed = grid_type(shape, {1.0, 1.0}, left_closed_status);
            grid_type raster_grid_vert_looped = grid_type(shape, {1.0, 1.0}, vert_looped_status);

            elev_type elevation
               {{0.5, 0.4, 3.0},
                {3.0, 0.1, 3.0},
                {3.0, 3.0, 3.0}};
        };

        class sinks_profile_grid: public ::testing::Test
        {
        protected:

            using grid_type = fs::profile_grid;
            using elev_type = xt::xtensor<double, 1>;

            fs::node_status ns_core = fs::node_status::core;
            fs::node_status ns_fixed = fs::node_status::fixed_value_boundary;

            grid_type profile_grid_closed = grid_type(4, 1.0, ns_fixed);
            grid_type profile_grid_half_open = grid_type(4, 1.0, {ns_fixed, ns_core});

            elev_type elevation {3.0, 0.1, 0.1, 2.0};
        };

    }
}

#endif // FASTSCAPELIB_TEST_SINKS_H
