#include "fastscape_benchmark.hpp"
#include "fastscapelib/bedrock_chanel.hpp"
#include "fastscapelib/flow_routing.hpp"
#include "fastscapelib/basin_graph.hpp"

#include "xtensor/xrandom.hpp"

#include "dbg_output.hpp"

void
fastscape_run(size_t nrows, size_t ncols, FastscapeFunctionType func)
{
    const double dx = 100.0;
    const double dy = 100.0;
    const double cell_area = dx * dy;
    const double uplift_rate = 2.e-3;
    const double k = 7e-4;
    const double n = 1.0;
    const double m = 0.4;
    const double dt = 10000.0;
    const double tolerance = 1e-3;

    const int nsteps = 30;

    std::array<size_t, 2> shape = { nrows, ncols };
    std::array<size_t, 1> shape1D = { nrows * ncols };

    xt::xtensor<double, 2> elevation = xt::random::rand(shape, 0.0, 1.0);

    xt::xtensor<bool, 2> active_nodes(elevation.shape());

    for (size_t i = 0; i < active_nodes.shape()[0]; ++i)
        for (size_t j = 0; j < active_nodes.shape()[1]; ++j)
            active_nodes(i, j) = i != 0 && j != 0 && i != active_nodes.shape()[0] - 1
                                 && j != active_nodes.shape()[1] - 1;

    double variation = std::min(1.0, 64.0 * 64.0 / double(nrows * ncols));

    xt::xtensor<double, 2> uplift = (1 - variation + variation * elevation);

    elevation = 1e-3 * elevation;


    xt::xtensor<index_t, 1> stack(shape1D);
    xt::xtensor<index_t, 1> receivers(shape1D);
    xt::xtensor<double, 1> dist2receivers(shape1D);
    xt::xtensor<double, 1> area(shape1D);

    xt::xtensor<double, 1> erosion(shape1D);

    // if needed, for computing the number of pits
    xt::xtensor<index_t, 1> pits(shape1D);
    xt::xtensor<index_t, 1> basins(shape1D);
    xt::xtensor<index_t, 1> ndonors(shape1D);
    xt::xtensor<index_t, 2> donors({ nrows * ncols, 8 });

    fs::BasinGraph<index_t, index_t, double> basin_graph;

    static int experiment = -1;
    ++experiment;
    std::stringstream ss;
    ss << "out/elevation-" << experiment << "-";

    std::stringstream ssd;
    ssd << "out/drain-" << experiment << "-";

    int first, last;


    for (int step_count = 0; step_count < nsteps; ++step_count)
    {
        if (nrows < 1025)
            dbg_out(ss.str(), step_count, elevation, shape);

        // apply uplift
        elevation += active_nodes * dt * uplift_rate * uplift;

        // std::cout << elevation << std::endl;

        if (true)
        {  // compute number of pits

            fs::compute_receivers_d8(receivers, dist2receivers, elevation, active_nodes, dx, dy);

            fs::compute_donors(ndonors, donors, receivers);
            fs::compute_stack(stack, ndonors, donors, receivers);

            basin_graph.compute_basins(basins, stack, receivers);
            auto outlets = xt::adapt(basin_graph.outlets(), shape1D);
            auto npits = fs::find_pits(pits, outlets, active_nodes, basin_graph.basin_count());

            std::cout << npits << ' ';

            if (step_count == 0)
                first = (int) npits;
            if (step_count == nsteps - 1)
                last = (int) npits;

            // std::cout << step_count << " " << npits << std::endl;

            // if( step_count && npits)
            //{
            //	std::cout << "err " << pits[0] << std::endl;
            // }

            /*if (step_count == 1)
            {
                for (int i = 35; i < 44; ++i)
                {
                    for (int j = 35; j < 44; ++j)
                    {
                        if (j + 44 * i == 1800)
                            std::cout << '[';
                        std::cout << elevation[j + 44 * i]-4;
                        if (j + 44 * i == 1800)
                            std::cout << ']';
                        std::cout << ' ';

                    }

                    std::cout << std::endl;
                }
            }*/
        }

        // compute stack
        func(stack, receivers, dist2receivers, elevation, active_nodes, dx, dy);


        // check if receivers are neighbors.
        if (false)
        {
            for (int i = 0; i < receivers.shape()[0]; ++i)
                if (active_nodes(i))
                {
                    int64_t rid = receivers(i);

                    auto x = i % ncols;
                    auto y = i / ncols;

                    auto rx = rid % ncols;
                    auto ry = rid / ncols;

                    if (x == rx && y == ry)
                    {
                        std::cout << i << "(" << x << "," << y << ") -> " << rid << std::endl;
                        std::abort();
                    }

                    if (x == rx)
                    {
                        if (dist2receivers(i) != dy)
                        {
                            std::cout << i << " " << dist2receivers(i) << ' ' << dx << ' ' << dy
                                      << std::endl;
                            std::abort();
                        }
                        assert(y == ry - 1 || y == ry + 1);
                    }
                    else if (y == ry)
                    {
                        assert(dist2receivers(i) == dx);
                        assert(x == rx - 1 || x == rx + 1);
                    }
                    else
                    {
                        assert(x == rx - 1 || x == rx + 1);
                        assert(y == ry - 1 || y == ry + 1);

                        if (dist2receivers(i) != std::sqrt(dx * dx + dy * dy))
                        {
                            std::cout << i << " " << dist2receivers(i) << ' '
                                      << std::sqrt(dx * dx + dy * dy) << std::endl;
                            std::abort();
                        }
                    }
                }
        }

        // update drainage area
        fs::compute_drainage_area(area, stack, receivers, dx, dy);
        if (nrows < 1025)
            dbg_out(ssd.str(), step_count, area, shape);

        //        double mx = std::numeric_limits<double>::lowest();
        //        double mn = std::numeric_limits<double>::max();

        //        for(int i = 0; i< area.shape()[0]; ++i)
        //        {
        //            mx = std::max(mx, area[i]);
        //            mn = std::min(mn, area[i]);
        //        }
        //        std::cout << mn << ' ' << mx << std::endl;

        fs::erode_spower(
            erosion, elevation, stack, receivers, dist2receivers, area, k, m, n, dt, tolerance);

        // apply erosion
        for (size_t k = 0; k < nrows * ncols; ++k)
            elevation(k) -= erosion(k);
    }

    std::cout << double(first) / double(nrows * ncols) << ' '
              << double(last) / double(nrows * ncols) << std::endl;
}
