#include <cstddef>
#include <array>
#include <stdexcept>

#include "gtest/gtest.h"
#include "xtensor/xtensor.hpp"
#include "xtensor/xoperation.hpp"
#include "xtensor/xview.hpp"

#include "fastscapelib/boundary.hpp"


namespace fs = fastscapelib;


TEST(boundary, create_node_status)
{
    std::array<std::size_t, 2> shape {5, 5};
    xt::xtensor<fs::NodeStatus, 2> node_status = fs::create_node_status(shape);

    EXPECT_TRUE(xt::all(xt::equal(node_status, fs::NodeStatus::CORE_NODE)));
}


TEST(boundary, set_node_status_grid_boundaries)
{
    auto node_status = xt::empty<fs::NodeStatus>({5, 5});

    {
        SCOPED_TRACE("test basic logic");

        fs::set_node_status_grid_boundaries(node_status,
                                            fs::NodeStatus::FIXED_VALUE_BOUNDARY,
                                            fs::NodeStatus::FIXED_VALUE_BOUNDARY,
                                            fs::NodeStatus::FIXED_VALUE_BOUNDARY,
                                            fs::NodeStatus::FIXED_VALUE_BOUNDARY);

        auto top_bottom_v = xt::view(node_status, xt::keep(0, -1), xt::all());

        EXPECT_TRUE(
            xt::all(xt::equal(top_bottom_v, fs::NodeStatus::FIXED_VALUE_BOUNDARY))
        );

        auto left_right_v = xt::view(node_status, xt::all(), xt::keep(0, -1));

        EXPECT_TRUE(
            xt::all(xt::equal(left_right_v, fs::NodeStatus::FIXED_VALUE_BOUNDARY))
            );
    }

    {
        SCOPED_TRACE("test exception non-symmetrical looped boundaries");

        EXPECT_THROW(
            fs::set_node_status_grid_boundaries(node_status,
                                                fs::NodeStatus::LOOPED_BOUNDARY,
                                                fs::NodeStatus::LOOPED_BOUNDARY,
                                                fs::NodeStatus::CORE_NODE,
                                                fs::NodeStatus::CORE_NODE),
            std::invalid_argument
            );

        EXPECT_THROW(
            fs::set_node_status_grid_boundaries(node_status,
                                                fs::NodeStatus::LOOPED_BOUNDARY,
                                                fs::NodeStatus::CORE_NODE,
                                                fs::NodeStatus::CORE_NODE,
                                                fs::NodeStatus::CORE_NODE),
            std::invalid_argument
            );

        EXPECT_THROW(
            fs::set_node_status_grid_boundaries(node_status,
                                                fs::NodeStatus::CORE_NODE,
                                                fs::NodeStatus::LOOPED_BOUNDARY,
                                                fs::NodeStatus::CORE_NODE,
                                                fs::NodeStatus::CORE_NODE),
            std::invalid_argument
            );
    }

    {
        SCOPED_TRACE("test solve conflicts at grid corners");

        fs::set_node_status_grid_boundaries(node_status,
                                            fs::NodeStatus::CORE_NODE,
                                            fs::NodeStatus::FIXED_VALUE_BOUNDARY,
                                            fs::NodeStatus::FIXED_GRADIENT_BOUNDARY,
                                            fs::NodeStatus::CORE_NODE);

        EXPECT_TRUE(node_status(0, 4) == fs::NodeStatus::FIXED_VALUE_BOUNDARY);
        EXPECT_TRUE(node_status(4, 4) == fs::NodeStatus::FIXED_VALUE_BOUNDARY);
        EXPECT_TRUE(node_status(4, 0) == fs::NodeStatus::FIXED_GRADIENT_BOUNDARY);
        EXPECT_TRUE(node_status(0, 0) == fs::NodeStatus::CORE_NODE);

        fs::set_node_status_grid_boundaries(node_status,
                                            fs::NodeStatus::CORE_NODE,
                                            fs::NodeStatus::LOOPED_BOUNDARY,
                                            fs::NodeStatus::FIXED_GRADIENT_BOUNDARY,
                                            fs::NodeStatus::LOOPED_BOUNDARY);

        EXPECT_TRUE(node_status(0, 4) == fs::NodeStatus::LOOPED_BOUNDARY);
        EXPECT_TRUE(node_status(4, 4) == fs::NodeStatus::FIXED_GRADIENT_BOUNDARY);
        EXPECT_TRUE(node_status(4, 0) == fs::NodeStatus::FIXED_GRADIENT_BOUNDARY);
        EXPECT_TRUE(node_status(0, 0) == fs::NodeStatus::LOOPED_BOUNDARY);
    }
}
