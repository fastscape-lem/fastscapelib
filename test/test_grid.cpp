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

        fs::edge_nodes_status es2 {fs::node_status::fixed_value_boundary,
                                   fs::node_status::fixed_gradient_boundary};

        EXPECT_EQ(es2.left, fs::node_status::fixed_value_boundary);
        EXPECT_EQ(es2.right, fs::node_status::fixed_gradient_boundary);

        EXPECT_THROW(fs::edge_nodes_status es3 {fs::node_status::core},
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

        fs::border_nodes_status bs2 {fs::node_status::fixed_value_boundary,
                                     fs::node_status::looped_boundary,
                                     fs::node_status::fixed_value_boundary,
                                     fs::node_status::looped_boundary};

        EXPECT_EQ(bs2.top, fs::node_status::fixed_value_boundary);
        EXPECT_EQ(bs2.right, fs::node_status::looped_boundary);
        EXPECT_EQ(bs2.bottom, fs::node_status::fixed_value_boundary);
        EXPECT_EQ(bs2.left, fs::node_status::looped_boundary);

        EXPECT_THROW(fs::border_nodes_status bs3 {fs::node_status::core},
                     std::invalid_argument);
    }
}

TEST(grid, profile_grid_constructor)
{
    auto es = fs::edge_nodes_status(fs::node_status::fixed_value_boundary);
    auto g = fs::profile_grid(10, 2.0, es);

    EXPECT_EQ(g.size(), size_t(10));
    EXPECT_EQ(g.spacing(), 2.0);
    EXPECT_EQ(g.status_at_nodes().size(), size_t(10));
    EXPECT_EQ(g.status_at_nodes()(0), fs::node_status::fixed_value_boundary);
    EXPECT_EQ(g.status_at_nodes()(5), fs::node_status::core);
    EXPECT_EQ(g.status_at_nodes()(9), fs::node_status::fixed_value_boundary);

    {
        SCOPED_TRACE("test status_at_nodes parameter");

        auto g2 = fs::profile_grid(10, 2.0, es,
                                   {{5, fs::node_status::fixed_value_boundary}});

        EXPECT_EQ(g2.status_at_nodes()(5), fs::node_status::fixed_value_boundary);
    }

    {
        SCOPED_TRACE("test invalid looped boundary status");

        fs::edge_nodes_status es2 {fs::node_status::fixed_value_boundary,
                                   fs::node_status::looped_boundary};

        EXPECT_THROW(fs::profile_grid(10, 2.0, es2), std::invalid_argument);
    }
}

TEST(grid, profile_grid_neighbors)
{
    auto es = fs::edge_nodes_status(fs::node_status::fixed_value_boundary);
    auto g = fs::profile_grid(10, 2.0, es);

    auto& n0 = g.neighbors(0);
    EXPECT_EQ(n0.size(), size_t(1));
    EXPECT_EQ(n0[0].idx, size_t(1));
    EXPECT_EQ(n0[0].distance, 2.0);
    EXPECT_EQ(n0[0].status, fs::node_status::core);

    auto& n9 = g.neighbors(9);
    EXPECT_EQ(n9.size(), size_t(1));
    EXPECT_EQ(n9[0].idx, size_t(8));

    auto& n5 = g.neighbors(5);
    EXPECT_EQ(n5.size(), size_t(2));
    EXPECT_EQ(n5[0].idx, size_t(4));
    EXPECT_EQ(n5[1].idx, size_t(6));

    {
        SCOPED_TRACE("test looped boundary neighbors");

        auto es2 = fs::edge_nodes_status(fs::node_status::looped_boundary);
        auto g2 = fs::profile_grid(10, 2.0, es2);

        auto& n20 = g2.neighbors(0);
        EXPECT_EQ(n20.size(), size_t(2));
        EXPECT_EQ(n20[0].idx, size_t(9));
        EXPECT_EQ(n20[1].idx, size_t(1));

        auto& n29 = g2.neighbors(9);
        EXPECT_EQ(n29.size(), size_t(2));
        EXPECT_EQ(n29[0].idx, size_t(8));
        EXPECT_EQ(n29[1].idx, size_t(0));
    }
}
