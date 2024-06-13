#include <iostream>

#include "xtensor/xarray.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xview.hpp"

#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/grid/raster_grid.hpp"
#include "fastscapelib/eroders/diffusion_adi.hpp"
#include "fastscapelib/eroders/spl.hpp"


namespace fs = fastscapelib;


int
main()
{
    // raster grid and boundary conditions
    fs::raster_boundary_status bs(fs::node_status::fixed_value);

    auto grid = fs::raster_grid::from_length({ 201, 301 }, { 5e4, 7.5e4 }, bs);

    // flow graph with single direction flow routing
    fs::flow_graph<fs::raster_grid> flow_graph(
        grid, { fs::single_flow_router(), fs::mst_sink_resolver() });

    // Setup eroders
    auto spl_eroder = fs::make_spl_eroder(flow_graph, 2e-4, 0.4, 1, 1e-5);
    auto diffusion_eroder = fs::diffusion_adi_eroder(grid, 0.01);

    // initial topographic surface elevation (flat surface + random perturbations)
    xt::xarray<double> init_elevation = xt::random::rand<double>(grid.shape());
    xt::xarray<double> elevation = init_elevation;

    // init drainage area and temp arrays
    xt::xarray<double> drainage_area(elevation.shape());
    xt::xarray<double> uplifted_elevation(elevation.shape());

    // uplift rate
    xt::xarray<double> uplift_rate(elevation.shape(), 1e-3);
    auto row_bounds = xt::view(uplift_rate, xt::keep(0, -1), xt::all());
    row_bounds = 0.0;
    auto col_bounds = xt::view(uplift_rate, xt::all(), xt::keep(0, -1));
    col_bounds = 0.0;

    //
    // run model
    //
    double dt = 2e4;
    int nsteps = 50;

    for (int step = 0; step < nsteps; step++)
    {
        // apply uplift
        uplifted_elevation = elevation + dt * uplift_rate;

        // flow routing
        flow_graph.update_routes(uplifted_elevation);

        // flow accumulation (drainage area)
        flow_graph.accumulate(drainage_area, 1.0);

        // apply channel erosion then hillslope diffusion
        auto spl_erosion = spl_eroder.erode(uplifted_elevation, drainage_area, dt);
        auto diff_erosion = diffusion_eroder.erode(uplifted_elevation - spl_erosion, dt);

        // update topography
        elevation = uplifted_elevation - spl_erosion - diff_erosion;
    }

    std::cout << "mean final elevation: " << xt::mean(elevation) << std::endl;
    std::cout << "max final elevation: " << xt::amax(elevation) << std::endl;

    return 0;
}
