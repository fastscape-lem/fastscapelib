#include "fastscape_benchmark.hpp"
#include "bedrock_chanel.hpp"

void fastscape_run(size_t nrows, size_t ncols, FastscapeFunctionType func)
{
    std::array<size_t, 2> shape = { nrows, ncols};

    xt::xtensor<double, 2> elevation = xt::random::rand(shape, 0.0, 1.0);

    xt::xtensor<bool, 2> active_nodes(elevation.shape());
    for(size_t i = 0; i<active_nodes.shape()[0]; ++i)
        for(size_t j = 0; j<active_nodes.shape()[1]; ++j)
            active_nodes(i,j) = i != 0 && j != 0
                    && i != active_nodes.shape()[0] -1
                    && j != active_nodes.shape()[1] -1;

    double uplift_tate = 1.0;
    double k = 1.0;
    double m = 1.0;
    double n = 1.0;
    double dt = 1.0;
    double tolerance = 1.0;


    // apply uplift
    elevation += active_nodes

    func(elevation, active_nodes); //->stack, receivers, dist2receivers

    //area
    compute_drainage_area(area, stack, receivers);

    erode_spower(erosion, elevation, stack, receivers, dist2receivers, area,
                 k, m, n, dt, tolerance);

    // apply erosion
}
