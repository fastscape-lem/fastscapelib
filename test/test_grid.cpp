#include "gtest/gtest.h"

#include "fastscapelib/grid.hpp"


namespace fs = fastscapelib;

TEST(grid, edge_nodes_status)
{
    {
        SCOPED_TRACE("single status constructor");

        auto es = fs::edge_nodes_status(fs::node_status::fixed_value_boundary);

        EXPECT_EQ(es.left, fs::node_status::fixed_value_boundary);
        EXPECT_EQ(es.right, fs::node_status::fixed_value_boundary);
    }

    {
        SCOPED_TRACE("initializer constructor");

        fs::edge_nodes_status es {fs::node_status::fixed_value_boundary,
                                  fs::node_status::fixed_gradient_boundary};

        EXPECT_EQ(es.left, fs::node_status::fixed_value_boundary);
        EXPECT_EQ(es.right, fs::node_status::fixed_gradient_boundary);

        EXPECT_THROW(fs::edge_nodes_status es2 {fs::node_status::core},
                     std::invalid_argument);
    }
}

TEST(grid, border_nodes_status)
{
    {
        SCOPED_TRACE("single status constructor");

        auto bs = fs::border_nodes_status(fs::node_status::fixed_value_boundary);

        EXPECT_EQ(bs.top, fs::node_status::fixed_value_boundary);
        EXPECT_EQ(bs.right, fs::node_status::fixed_value_boundary);
        EXPECT_EQ(bs.bottom, fs::node_status::fixed_value_boundary);
        EXPECT_EQ(bs.left, fs::node_status::fixed_value_boundary);
    }

    {
        SCOPED_TRACE("initializer constructor");

        fs::border_nodes_status bs {fs::node_status::fixed_value_boundary,
                                    fs::node_status::looped_boundary,
                                    fs::node_status::fixed_value_boundary,
                                    fs::node_status::looped_boundary};

        EXPECT_EQ(bs.top, fs::node_status::fixed_value_boundary);
        EXPECT_EQ(bs.right, fs::node_status::looped_boundary);
        EXPECT_EQ(bs.bottom, fs::node_status::fixed_value_boundary);
        EXPECT_EQ(bs.left, fs::node_status::looped_boundary);

        EXPECT_THROW(fs::border_nodes_status bs2 {fs::node_status::core},
                     std::invalid_argument);
    }
}
