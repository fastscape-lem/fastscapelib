#include "dbg_output.hpp"
#include "fastscapelib/basin_graph.hpp"
#include "fastscapelib/flow_routing.hpp"
#include "fastscapelib/bedrock_chanel.hpp"

#include "examples.hpp"


template <fs::ConnectType connect, class Elev_T, class Active_T>
void fastscape(Elev_T& elevation, const Active_T& active_nodes)
{
	auto shape = elevation.shape();
	int nrows = (int)shape[0];
	int ncols = (int)shape[1];
	std::array<size_t, 1> shape1D = { (size_t)nrows * (size_t)ncols };

	double dx = 100.0, dy = 100.0;

	fastscapelib::BasinGraph<index_t, index_t, double> bg;
	xt::xtensor<index_t, 1> basins(shape1D);
	xt::xtensor<index_t, 1> stack(shape1D);
	xt::xtensor<index_t, 1> receivers(shape1D);
	xt::xtensor<double, 1> dist2receivers(shape1D);
	xt::xtensor<index_t, 1> ndonors(shape1D);
	xt::xtensor<index_t, 2> donors({ shape1D[0], 8 });
	xt::xtensor<double, 1> area(shape1D);
	xt::xtensor<double, 1> erosion(shape1D);

	int num_iter = 10;

	for (int s = 0; s < num_iter; ++s)
	{

		fastscapelib::compute_receivers_d8(receivers, dist2receivers,
			elevation, active_nodes,
			dx, dy);

		fastscapelib::compute_donors(ndonors, donors, receivers);
		fastscapelib::compute_stack(stack, ndonors, donors, receivers);

		fastscapelib::correct_flowrouting<fs::BasinAlgo::Boruvka, fs::ConnectType::Sloped>(bg, basins, receivers, dist2receivers,
			ndonors, donors, stack, active_nodes, elevation, dx, dy);

		if (s != num_iter - 1)
		{
			area = xt::ones<index_t>({ nrows*ncols }) * dx*dy;
			fs::compute_drainage_area(area, stack, receivers);
			fs::erode_spower(erosion, elevation, stack, receivers, dist2receivers, area,
				7.0e-4f, 0.4f, 1.0f, 10.0f, 1.0e-4f);

			for (size_t k = 0; k < nrows*ncols; ++k)
				elevation(k) = elevation(k) - erosion(k);
		}
		

	}
}

void example_jail()
{

	int nrows = 500, ncols = 500;
	
	std::array<size_t, 2> shape = { (size_t)nrows, (size_t)ncols };

	xt::xtensor<double, 2> elevation(shape);
	xt::xtensor<double, 2> drain(shape);

	for (int y = 0; y<nrows; ++y)
		for (int x = 0; x < ncols; ++x)
			elevation(y, x) = 0.5 * double(std::abs(x - (ncols / 2)) + std::abs(y - (nrows / 2)));

	elevation(0, ncols / 2) = 0.0;
	
	xt::xtensor<bool, 2> active_nodes = xt::ones<bool>(elevation.shape());
	active_nodes(0, ncols / 2) = false;

	xt::xtensor<double, 2> elev_simple = elevation;
	fastscape<fs::ConnectType::Simple>(elev_simple, active_nodes);
	xt::xtensor<double, 2> elev_carved = elevation;
	fastscape<fs::ConnectType::Carved>(elev_carved, active_nodes);
	xt::xtensor<double, 2> elev_sloped = elevation;
	fastscape<fs::ConnectType::Sloped>(elev_sloped, active_nodes);


	dbg_out("results/jail/elevation-", 0, elevation, shape);
	dbg_out("results/jail/simple-", 0, elev_simple, shape);
	dbg_out("results/jail/carved-", 0, elev_carved, shape);
	dbg_out("results/jail/sloped-", 0, elev_sloped, shape);

}