#include <memory>

#include "xtensor/xarray.hpp"
#include "xtensor/xbuilder.hpp"

#include "gtest/gtest.h"

#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/grid/raster_grid.hpp"

#include "flow_operator_fixtures.hpp"


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace testing
    {

        TEST(flow_operator, default_member_values)
        {
            default_operator op;

            EXPECT_EQ(op.name(), "default_operator");
            EXPECT_EQ(op.graph_updated, false);
            EXPECT_EQ(op.elevation_updated, false);
            EXPECT_EQ(op.in_flowdir, fs::flow_direction::undefined);
            EXPECT_EQ(op.out_flowdir, fs::flow_direction::undefined);
        }

        TEST(flow_operator, return_same_elevation_ref)
        {
            fs::raster_boundary_status bs{ fs::node_status::fixed_value_boundary };
            auto grid = fs::raster_grid({ 5, 5 }, { 1.0, 1.0 }, bs);
            xt::xarray<double> elevation = xt::zeros<double>(grid.shape());

            test_operator op;
            auto graph = fs::flow_graph<fs::raster_grid>(grid, { fake_operator() });

            const auto& actual = graph.update_routes(elevation);

            // when elevation is not updated by any of the operators, update_routes
            // should return a const reference of the input elevation.
            ASSERT_TRUE(&actual == &elevation);
        }

        TEST(flow_operator, apply)
        {
            fs::raster_boundary_status bs{ fs::node_status::fixed_value_boundary };
            auto grid = fs::raster_grid({ 5, 5 }, { 1.0, 1.0 }, bs);
            xt::xarray<double> elevation = xt::zeros<double>(grid.shape());

            auto op = std::make_shared<test_operator>();
            auto graph = fs::flow_graph<fs::raster_grid>(grid, { op, fake_operator() });

            {
                SCOPED_TRACE("test apply 1st pass");

                const auto& actual = graph.update_routes(elevation);
                xt::xarray<double> expected = elevation + op->m_offset;

                EXPECT_EQ(actual, expected);
            }

            {
                SCOPED_TRACE("test apply 2nd pass (update operator params)");

                op->m_offset = 2.0;
                const auto& actual = graph.update_routes(elevation);
                xt::xarray<double> expected = elevation + op->m_offset;

                EXPECT_EQ(actual, expected);
            }
        }
    }
}
