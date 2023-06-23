#include <iostream>

#include "xtensor/xtensor.hpp"
#include "xtensor/xmath.hpp"

#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/grid/profile_grid.hpp"
#include "fastscapelib/eroders/spl.hpp"


namespace fs = fastscapelib;


int
main()
{
    // profile grid and boundary conditions
    using ns = fs::node_status;
    fs::profile_boundary_status bs{ ns::fixed_value, ns::core };

    fs::profile_grid grid = fs::profile_grid(101, 300.0, bs);

    // grid x-coordinate values
    double x0 = 300.;
    double length = (grid.size() - 1.0) * grid.spacing();
    xt::xtensor<double, 1> x = xt::linspace<double>(length + x0, x0, grid.size());

    // flow graph with single direction flow routing
    fs::flow_graph<fs::profile_grid> flow_graph(grid, { fs::single_flow_router() });

    // compute drainage area using Hack's law
    double hack_coef = 6.69;
    double hack_exp = 1.67;
    xt::xtensor<double, 1> drainage_area = hack_coef * xt::pow(x, hack_exp);

    // Setup the SPL eroder
    auto spl_eroder = fs::make_spl_eroder(flow_graph, 1e-4, 0.5, 1, 1e-5);

    // initial river profile elevation (gentle linear slope)
    xt::xtensor<double, 1> init_elevation = (length + x0 - x) * 1e-4;
    xt::xtensor<double, 1> elevation = init_elevation;

    // uplift rate
    xt::xtensor<double, 1> uplift_rate(x.shape(), 1e-3);
    uplift_rate[0] = 0.0;

    //
    // run model
    //
    double dt = 1e2;
    int nsteps = 4000;

    // compute flow paths (needed only once in this example)
    flow_graph.update_routes(elevation);

    for (int step = 0; step < nsteps; step++)
    {
        // apply uplift
        xt::xtensor<double, 1> uplifted_elevation = elevation + dt * uplift_rate;

        // abrupt change in erodibility
        if (step == 3000)
        {
            spl_eroder.set_k_coef(spl_eroder.k_coef() / 4.0);
        }

        // apply channel erosion
        auto spl_erosion = spl_eroder.erode(elevation, drainage_area, dt);

        // update topography
        elevation = uplifted_elevation - spl_erosion;
    }

    std::cout << "mean final elevation: " << xt::mean(elevation) << std::endl;
    std::cout << "max final elevation: " << xt::amax(elevation) << std::endl;

    return 0;
}
