#include "fastscape_benchmark.hpp"
#include "fastscapelib\bedrock_chanel.hpp"
#include "fastscapelib\flow_routing.hpp"

#include "xtensor/xrandom.hpp"

void fastscape_run(size_t nrows, size_t ncols, FastscapeFunctionType func)
{
	const double dx = 100.0;
	const double dy = 100.0;
	const double cell_area = dx*dy;
	const double uplift_rate = 2.e-3;
	const double k = 7e-4;
	const double n = 1.0;
	const double m = 0.4;
	const double dt = 1000.0;
	const double tolerance = 1e-3;

	const int nsteps = 50;

	std::array<size_t, 2> shape = { nrows, ncols };
	std::array<size_t, 1> shape1D = { nrows * ncols };

    xt::xtensor<double, 2> elevation = xt::random::rand(shape, 0.0, 1.0);

    xt::xtensor<bool, 2> active_nodes(elevation.shape());

    for(size_t i = 0; i<active_nodes.shape()[0]; ++i)
        for(size_t j = 0; j<active_nodes.shape()[1]; ++j)
            active_nodes(i,j) = i != 0 && j != 0
                    && i != active_nodes.shape()[0] -1
                    && j != active_nodes.shape()[1] -1;

	xt::xtensor<index_t, 1> stack (shape1D);
	xt::xtensor<index_t, 1> receivers(shape1D);
	xt::xtensor<double, 1> dist2receivers(shape1D);
	xt::xtensor<double, 1> area(shape1D);

	xt::xtensor<double, 1> erosion(shape1D);


	for (int i = 0; i < nsteps; ++i)
	{
		// apply uplift
		elevation += active_nodes * dt * uplift_rate;

		// compute stack
		func(stack, receivers, dist2receivers, elevation, active_nodes, dx, dy);

		// update drainage area
		area = xt::ones<index_t>({ nrows*ncols }) * cell_area;
		fs::compute_drainage_area(area, stack, receivers);

		fs::erode_spower(erosion, elevation, stack, receivers, dist2receivers, area,
			k, m, n, dt, tolerance);

		// apply erosion
		for(size_t i = 0; i< nrows*ncols; ++i)
			elevation (i) -= erosion(i);
	}
}
