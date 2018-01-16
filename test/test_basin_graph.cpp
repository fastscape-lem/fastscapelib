#include <cmath>

#include "gtest/gtest.h"
#include "xtensor/xtensor.hpp"
#include "xtensor/xio.hpp"

#include "fastscapelib/utils.hpp"
#include "fastscapelib/flow_routing.hpp"
#include "fastscapelib/basin_graph.hpp"


class BasinGraph_Test
{
public:
    BasinGraph_Test()
    {

        elevation = xt::xtensor<double, 2>{{2.0, 3.0, 3.0, 2.0, 5.0},
                                           {8.0, 0.0, 1.0, 1.0, 4.0},
                                           {2.0, 4.0, 3.0, 2.0, 6.0}};
        active_nodes =  xt::xtensor<bool, 2>{{false, false, false, false, false},
                                             {false, true,  true,  true, false},
                                             {false, false, false, false, false}};


        receivers = xt::ones<index_t>({15}) * -1;
        dist2receivers = xt::ones<double>({15}) * -1.;
        fs::compute_receivers_d8(receivers, dist2receivers,
                                 elevation, active_nodes,
                                 1., 1.);

        ndonors = xt::ones<index_t>({15}) * -1;
        donors = xt::ones<index_t>({15, 8}) * -1;

        fs::compute_donors(ndonors, donors, receivers);

        stack = xt::ones<index_t>({15}) * -1;
        fs::compute_stack(stack, ndonors, donors, receivers);

        basins = xt::ones<index_t>({15}) * -1;
    }


    void test_basins()
    {
        /////////// Answers ///////////
        xt::xtensor<index_t, 1> expected_basins
        {0,  1,  2,  3,  4,
            5,  6,  6,  7,  8,
            9, 10, 11, 12, 13};
        xt::xtensor<index_t, 1> expected_outlets
        {0,  1,  2,  3,  4,
            5,  6,  8,  9,
            10, 11, 12, 13, 14};


        /////////// run test case ///////////
        xt::xtensor<index_t, 1> basins = xt::ones<index_t>({15}) * -1;

        basin_graph.compute_basins(basins, stack, receivers);
        auto& outlets = basin_graph.outlets();

        EXPECT_EQ(outlets.size(), expected_outlets.shape()[0]);
        EXPECT_TRUE(xt::all(xt::equal(basins, expected_basins)));

        EXPECT_TRUE(xt::all(xt::equal(xt::adapt(outlets, expected_outlets.shape()), expected_outlets)));
    }

    void test_connect()
    {
        double minf = std::numeric_limits<double>::lowest();

        std::vector<fs::BasinGraph<index_t, index_t>::Link_T> expected_links =
        { {{0,1},{-1,-1},minf},
          {{0,2},{-1,-1},minf},
          {{0,3},{-1,-1},minf},
          {{0,4},{-1,-1},minf},
          {{0,5},{-1,-1},minf},
          {{6,1},{6,1},3},
          {{6,5},{6,5},8},
          {{6,10},{6,11},4},
          {{6,2},{7,2},3},
          {{6,7},{7,8},1},
          {{6,11},{7,12},3},
          {{7,3},{8,3},2},
          {{7,8},{8,9},4},
          {{7,12},{8,13},2},
          {{0,8},{-1,-1},minf},
          {{0,9},{-1,-1},minf},
          {{0,10},{-1,-1},minf},
          {{0,11},{-1,-1},minf},
          {{0,12},{-1,-1},minf},
          {{0,13},{-1,-1},minf}};


        /////////// run test case ///////////

        basin_graph.compute_basins(basins, stack, receivers);
        basin_graph.connect_basins(basins, receivers, stack, active_nodes, elevation);

        EXPECT_TRUE(is_links_eq(expected_links));
    }

    template<class Link_T>
    bool is_links_eq(std::vector<Link_T>& oth)
    {
        if (oth.size() != basin_graph._links.size())
            return false;
        for (size_t i = 0; i< basin_graph._links.size(); ++i)
            if (!(basin_graph._links[i] == oth[i]))
            {
                std::cout << '(' << basin_graph._links[i].basins[0]<<'-'<<basin_graph._links[i].basins[1] << "),(" << basin_graph._links[i].nodes[0]<<'-'<<basin_graph._links[i].nodes[1]  << ") w=" << basin_graph._links[i].weight << std::endl;
                std::cout << '(' << oth[i].basins[0]<<'-'<<oth[i].basins[1] << "),(" << oth[i].nodes[0]<<'-'<<oth[i].nodes[1]  << ") w=" << oth[i].weight << std::endl;
                return false;
            }
        return true;
    }
    void print_links()
    {
        for (auto& l : basin_graph._links)
            std::cout << '(' << l.basins[0]<<'-'<<l.basins[1] << "),(" << l.nodes[0]<<'-'<<l.nodes[1]  << ") w=" << l.weight << std::endl;
    }

    fs::BasinGraph<index_t, index_t> basin_graph;

    xt::xtensor<double, 2> elevation;
    xt::xtensor<bool, 2> active_nodes;

    xt::xtensor<index_t, 1> receivers;
    xt::xtensor<double, 1> dist2receivers;
    xt::xtensor<index_t, 1> ndonors;
    xt::xtensor<index_t, 2> donors;
    xt::xtensor<index_t, 1> stack;

    xt::xtensor<index_t, 1> basins;

};

TEST(basin_graph, compute_basins_1)
{
    fs::BasinGraph<size_t, size_t> basin_graph;

    // Example in Braun and Willet, 2013 as a test case.
    // TODO: test with multiple basins in node network.
    xt::xtensor<index_t, 1> receivers {1, 4, 1, 6, 4, 4, 5, 4, 6, 7};
    xt::xtensor<index_t, 1> stack {4, 1, 0, 2, 5, 6, 3, 8, 7, 9};
    xt::xtensor<index_t, 1> basins = xt::ones<index_t>({10}) * -1;

    xt::xtensor<index_t, 1> expected_basins
    {0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
    xt::xtensor<index_t, 1> expected_outlets
    {4, -1, -1, -1, -1, -1, -1, -1, -1, -1};

    basin_graph.compute_basins(basins, stack, receivers);
    auto& outlets = basin_graph.outlets();


    EXPECT_EQ(outlets.size(), 1);
    EXPECT_EQ(outlets[0], 4);
    EXPECT_TRUE(xt::all(xt::equal(basins, expected_basins)));
}

TEST(basin_graph, compute_basins_2)
{

    BasinGraph_Test bgt;

    bgt.test_basins();

}

TEST(basin_graph, connect_basins)
{
    BasinGraph_Test bgt;

    bgt.test_connect();
}
