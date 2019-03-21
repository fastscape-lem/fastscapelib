#include <array>

#include "gtest/gtest.h"

#include "fastscapelib/grid.hpp"


namespace fs = fastscapelib;


TEST(grid, boundary_status_constructors)
{
    {
        SCOPED_TRACE("single status");

        auto bs1 = fs::boundary_status(fs::node_status::fixed_value_boundary);

        EXPECT_EQ(bs1.left, fs::node_status::fixed_value_boundary);
        EXPECT_EQ(bs1.right, fs::node_status::fixed_value_boundary);
        EXPECT_EQ(bs1.top, fs::node_status::fixed_value_boundary);
        EXPECT_EQ(bs1.bottom, fs::node_status::fixed_value_boundary);
    }

    {
        SCOPED_TRACE("2-length initializer list");

        fs::boundary_status bs2 {fs::node_status::fixed_value_boundary,
                                 fs::node_status::fixed_gradient_boundary};

        EXPECT_EQ(bs2.left, fs::node_status::fixed_value_boundary);
        EXPECT_EQ(bs2.right, fs::node_status::fixed_gradient_boundary);
    }

    {
        SCOPED_TRACE("4-length initializer list");

        fs::boundary_status bs4 = {fs::node_status::fixed_value_boundary,
                                   fs::node_status::fixed_gradient_boundary,
                                   fs::node_status::looped_boundary,
                                   fs::node_status::looped_boundary};

        EXPECT_EQ(bs4.left, fs::node_status::fixed_value_boundary);
        EXPECT_EQ(bs4.right, fs::node_status::fixed_gradient_boundary);
        EXPECT_EQ(bs4.top, fs::node_status::looped_boundary);
        EXPECT_EQ(bs4.bottom, fs::node_status::looped_boundary);
    }

    {
        SCOPED_TRACE("2-length array");

        std::array<fs::node_status, 2> a2 = {fs::node_status::fixed_value_boundary,
                                             fs::node_status::fixed_gradient_boundary};
        fs::boundary_status bsa2(a2);

        EXPECT_EQ(bsa2.left, fs::node_status::fixed_value_boundary);
        EXPECT_EQ(bsa2.right, fs::node_status::fixed_gradient_boundary);
    }

    {
        SCOPED_TRACE("4-length array");

        std::array<fs::node_status, 4> a4 = {fs::node_status::fixed_value_boundary,
                                             fs::node_status::fixed_gradient_boundary,
                                             fs::node_status::looped_boundary,
                                             fs::node_status::looped_boundary};
        fs::boundary_status bsa4(a4);

        EXPECT_EQ(bsa4.left, fs::node_status::fixed_value_boundary);
        EXPECT_EQ(bsa4.right, fs::node_status::fixed_gradient_boundary);
        EXPECT_EQ(bsa4.top, fs::node_status::looped_boundary);
        EXPECT_EQ(bsa4.bottom, fs::node_status::looped_boundary);
    }
}

TEST(grid, profile_grid_constructor)
{
    auto g = fs::profile_grid(10, 2.0,
                              fs::boundary_status(fs::node_status::fixed_value_boundary));

    EXPECT_EQ(g.size(), size_t(10));
    EXPECT_EQ(g.spacing(), 2.0);
}

TEST(grid, profile_grid_status_at_nodes)
{
    auto g = fs::profile_grid(10, 2.0,
                              {fs::node_status::fixed_value_boundary,
                               fs::node_status::fixed_value_boundary});

    EXPECT_EQ(g.status_at_nodes().size(), size_t(10));
    EXPECT_EQ(g.status_at_nodes()(0), fs::node_status::fixed_value_boundary);
    EXPECT_EQ(g.status_at_nodes()(5), fs::node_status::core);
    EXPECT_EQ(g.status_at_nodes()(9), fs::node_status::fixed_value_boundary);

    {
        SCOPED_TRACE("test status_at_nodes grid constructor parameter");

        auto g2 = fs::profile_grid(
            10, 2.0,
            fs::boundary_status(fs::node_status::fixed_value_boundary),
            {{5, fs::node_status::fixed_value_boundary}});

        EXPECT_EQ(g2.status_at_nodes()(5), fs::node_status::fixed_value_boundary);
    }

    {
        SCOPED_TRACE("test invalid looped boundary status");

        fs::boundary_status bs {fs::node_status::fixed_value_boundary,
                                fs::node_status::looped_boundary};

        EXPECT_THROW(fs::profile_grid(10, 2.0, bs), std::invalid_argument);
    }
}

TEST(grid, profile_grid_neighbors)
{
    auto g = fs::profile_grid(10, 2.0,
                              fs::boundary_status(fs::node_status::fixed_value_boundary));

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

        auto g2 = fs::profile_grid(
            10, 2.0, fs::boundary_status(fs::node_status::looped_boundary));

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
