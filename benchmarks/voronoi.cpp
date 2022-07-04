#include <iostream>

#include "fastscapelib/basin_graph.hpp"
#include "fastscapelib/flow_routing.hpp"
#include "fastscapelib/bedrock_chanel.hpp"
#include "xtensor/xrandom.hpp"

#include "examples.hpp"
#include "dbg_output.hpp"

void
example_vornoi()
{
    int nrows = 500, ncols = 500;

    int vspace = 100;

    std::array<size_t, 2> shape = { (size_t) nrows, (size_t) ncols };
    std::array<size_t, 1> shape1D = { (size_t) nrows * (size_t) ncols };

    std::array<size_t, 2> vshape = { size_t(nrows / vspace), size_t(ncols / vspace) };

    xt::xtensor<double, 2> voronoi_x = xt::random::rand(shape, 0.1, 0.9);
    xt::xtensor<double, 2> voronoi_y = xt::random::rand(shape, 0.1, 0.9);

    xt::xtensor<double, 2> elevation(shape);
    xt::xtensor<double, 2> drain(shape);


    for (int y = 0; y < ncols; ++y)
        for (int x = 0; x < nrows; ++x)
        {
            double d = std::numeric_limits<double>::max();

            int vx = x / vspace;
            int vy = y / vspace;
            for (int py = std::max(0, vy - 2); py <= std::min((int) vshape[0] - 1, vy + 2); ++py)
                for (int px = std::max(0, vx - 2); px <= std::min((int) vshape[1] - 1, vx + 2);
                     ++px)
                {
                    double dx = double(vspace) * ((double) px + voronoi_x(py, px)) - (double) x;
                    double dy = double(vspace) * ((double) py + voronoi_y(py, px)) - (double) y;

                    double td = dx * dx + dy * dy;
                    if (td < d)
                    {
                        d = td;
                    }
                }

            elevation(y, x) = std::sqrt(d);
            drain(y, x) = std::sqrt(d);
        }


    fastscapelib::BasinGraph<index_t, index_t, double> bg;
    xt::xtensor<index_t, 1> basins(shape1D);
    xt::xtensor<index_t, 1> stack(shape1D);
    xt::xtensor<index_t, 1> receivers(shape1D);
    xt::xtensor<double, 1> dist2receivers(shape1D);
    xt::xtensor<index_t, 1> ndonors(shape1D);
    xt::xtensor<index_t, 2> donors({ (size_t) (nrows * ncols), 8 });

    xt::xtensor<bool, 2> active_nodes(elevation.shape());

    for (size_t i = 0; i < active_nodes.shape()[0]; ++i)
        for (size_t j = 0; j < active_nodes.shape()[1]; ++j)
            active_nodes(i, j) = i != 0 && j != 0 && i != active_nodes.shape()[0] - 1
                                 && j != active_nodes.shape()[1] - 1;

    double dx = 100.0, dy = 100.0;
    fastscapelib::compute_receivers_d8(receivers, dist2receivers, elevation, active_nodes, dx, dy);

    fastscapelib::compute_donors(ndonors, donors, receivers);
    fastscapelib::compute_stack(stack, ndonors, donors, receivers);

    fastscapelib::correct_flowrouting<fs::BasinAlgo::Boruvka, fs::ConnectType::Carved>(
        bg,
        basins,
        receivers,
        dist2receivers,
        ndonors,
        donors,
        stack,
        active_nodes,
        elevation,
        dx,
        dy);


    xt::xtensor<double, 1> dbasins(shape1D);
    xt::xtensor<double, 1> links(shape1D);
    xt::xtensor<double, 1> tree(shape1D);

    auto npits = bg.outlets().size();

    for (size_t id = 0; id < shape1D[0]; ++id)
    {
        dbasins[id] = (double) basins[id];
    }

    tree = xt::zeros<double>(shape1D);
    links = xt::zeros<double>(shape1D);

    for (auto& l : bg.getLinks())
    {
        if (l.nodes[0] >= 0)
            links[l.nodes[0]] = 1.0;

        if (l.nodes[1] >= 0)
            links[l.nodes[1]] = 1.0;
    }

    for (auto& t : bg.getTreeIndices())
    {
        auto& l = bg.getLinks()[t];
        if (l.nodes[0] >= 0)
            tree[l.nodes[0]] = 1.0;
        if (l.nodes[1] >= 0)
            tree[l.nodes[1]] = 2.0;
    }

    dbg_out("results/voronoi/elevation-", 0, elevation, shape);
    dbg_out("results/voronoi/basins-", 0, dbasins, shape);
    dbg_out("results/voronoi/links-", 0, links, shape);
    dbg_out("results/voronoi/tree-", 0, tree, shape);


    xt::xtensor<double, 1> area(shape1D);
    xt::xtensor<double, 1> erosion(shape1D);

    for (int s = 0; s < 10; ++s)
    {
        fs::compute_drainage_area(area, stack, receivers, dx, dy);
        fs::erode_spower(erosion,
                         elevation,
                         stack,
                         receivers,
                         dist2receivers,
                         area,
                         7.0e-4,
                         0.4,
                         1.0,
                         10.0,
                         1.0e-4);

        for (size_t k = 0; k < nrows * ncols; ++k)
            elevation(k) = elevation(k) - erosion(k);


        fastscapelib::compute_receivers_d8(
            receivers, dist2receivers, elevation, active_nodes, dx, dy);

        fastscapelib::compute_donors(ndonors, donors, receivers);
        fastscapelib::compute_stack(stack, ndonors, donors, receivers);

        fastscapelib::correct_flowrouting<fs::BasinAlgo::Boruvka, fs::ConnectType::Sloped>(
            bg,
            basins,
            receivers,
            dist2receivers,
            ndonors,
            donors,
            stack,
            active_nodes,
            elevation,
            dx,
            dy);
    }

    xt::xtensor<double, 1> water(shape1D);
    fs::fill_sinks_flat(water, elevation, stack, receivers);

    dbg_out("results/voronoi/eroded-", 0, elevation, shape);
    dbg_out("results/voronoi/water-", 0, water, shape);
}
