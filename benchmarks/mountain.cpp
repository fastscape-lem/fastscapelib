#include "dbg_output.hpp"
#include "fastscapelib/basin_graph.hpp"
#include "fastscapelib/flow_routing.hpp"
#include "fastscapelib/bedrock_chanel.hpp"

#include <xtensor/xrandom.hpp>

#include "examples.hpp"

template <class Elev_T, class Active_T>
void fastscape(Elev_T& elevation, const Active_T& active_nodes, double dx, double dt, int num_iter)
{
	auto shape = elevation.shape();
	int nrows = (int)shape[0];
	int ncols = (int)shape[1];
	std::array<size_t, 1> shape1D = { (size_t)nrows * (size_t)ncols };

	double dy = dx;

	fastscapelib::BasinGraph<index_t, index_t, double> bg;
	xt::xtensor<index_t, 1> basins(shape1D);
	xt::xtensor<index_t, 1> stack(shape1D);
	xt::xtensor<index_t, 1> receivers(shape1D);
	xt::xtensor<double, 1> dist2receivers(shape1D);
	xt::xtensor<index_t, 1> ndonors(shape1D);
	xt::xtensor<index_t, 2> donors({ shape1D[0], 8 });
	xt::xtensor<double, 1> area(shape1D);
	xt::xtensor<double, 1> erosion(shape1D);


	for (int s = 0; s < num_iter; ++s)
	{

		// uplift
		elevation += xt::cast<double>(active_nodes) * 5e-3 *dt;

		fastscapelib::compute_receivers_d8(receivers, dist2receivers,
			elevation, active_nodes,
			dx, dy);

		fastscapelib::compute_donors(ndonors, donors, receivers);
		fastscapelib::compute_stack(stack, ndonors, donors, receivers);

		fastscapelib::correct_flowrouting<fs::BasinAlgo::Boruvka, fs::ConnectType::Carved>(bg, basins, receivers, dist2receivers,
			ndonors, donors, stack, active_nodes, elevation, dx, dy);

		area = xt::ones<index_t>({ nrows*ncols }) * dx*dy;
		fs::compute_drainage_area(area, stack, receivers);
		fs::erode_spower(erosion, elevation, stack, receivers, dist2receivers, area,
			7.0e-4, 0.4, 1.0, dt, 1.0e-4);

		for (size_t k = 0; k < nrows*ncols; ++k)
			elevation(k) = elevation(k) - erosion(k);
	}
}

template <class T>
typename std::decay_t<T>::value_type fetch(const T& a, double r, double c)
{

	using F = std::decay_t<T>::value_type;

	auto shape = a.shape();

	// r \in [0, nr-1]

	int ir = std::min((int) r, (int)shape[0] - 2);
	int ic = std::min((int) c, (int)shape[1] - 2);

	double sr = r - double(ir);
	double sc = c - double(ic);

	F a00 = a(ir, ic);
	F a01 = a(ir, ic+1);
	F a10 = a(ir+1, ic);
	F a11 = a(ir+1, ic+1);

	return a00 * (1.0 - sr) * (1.0 - sc) + a01 * (1.0 - sr) * sc + a10 * sr * (1.0 - sc) + a11 * sc*sr;
}

template <class T>
void sample(const T& in, T& out)
{
	for (int r = 0; r < out.shape()[0]; ++r)
		for (int c = 0; c < out.shape()[1]; ++c)
		{
			double dr = double(r) / double(out.shape()[0] - 1) * double(in.shape()[0] - 1);
			double dc = double(c) / double(out.shape()[1] - 1) * double(in.shape()[1] - 1);
			out(r, c) = fetch(in, dr, dc);
		}
}

void generate_mountain()
{
	int n = 32;
	double dx = 64 * 10;
	double dt = 1000000;


	xt::xtensor<double, 2> elevation = xt::random::rand({ (size_t)n, (size_t)n }, 0.0, 1e-3);

	while (n <= 4096)
	//while (n <= 0)
	{

		xt::xtensor<bool, 2> active_nodes(elevation.shape());

		for (size_t i = 0; i<active_nodes.shape()[0]; ++i)
			for (size_t j = 0; j<active_nodes.shape()[1]; ++j)
				active_nodes(i, j) = i != 0 && j != 0
				&& i != active_nodes.shape()[0] - 1
				&& j != active_nodes.shape()[1] - 1;

		int num_iter = n == 32 ? 100 : 4;

		fastscape(elevation, active_nodes, dx, dt, num_iter);

		dbg_out("results/mountain/bedrock-", n, elevation, elevation.shape());

		xt::xtensor<double, 2> tmp({ (size_t)(2 * n), (size_t)(2 * n) }, -1);
		
		// linear upsampling
		sample(elevation, tmp);

		elevation = tmp;
		n *= 2;
		dx *= 0.5;
		dt *= 0.5;
	}

}

template <class S>
std::vector<size_t> holes(const S& shape)
{
	std::vector<size_t> ids;
	ids.reserve((shape[0] - 1) / 2 * ((shape[1] - 1) / 2));
	for (size_t r = 1; r < shape[0] - 1; r += 2)
		for (size_t c = 1; c < shape[1] - 1; c += 2)
			ids.push_back(c + r * shape[1]);

	assert(ids.size() == (shape[0] - 1) / 2 * ((shape[1] - 1) / 2));

	std::mt19937 gen = std::mt19937(123); // seed = 123
	std::shuffle(ids.begin(), ids.end(), gen);

	return ids;
}

template <class T>
void dig_hole(T& a, size_t hid)
{
	size_t r = hid / a.shape()[1];
	size_t c = hid % a.shape()[1];

	double h = std::numeric_limits<double>::max();
	for (size_t rr = r-1; rr <=r+1; ++rr)
		for (size_t cc = c - 1; cc <= c + 1; ++cc)
			h = std::min(h, a(rr, cc));

	a(r, c) = h - 1e-5;
}



void example_mountain()
{
	xt::xtensor<double, 1> h_prop = { 0, .015625, .03125, .0625, .125, .25, .5, 1.0 };
	xt::xtensor<int, 1> m_size = { 32, 64, 128, 256, 512, 1024, 2048, 4096 };

	xt::xtensor<double, 2> results({ m_size.size(), h_prop.size() });

	xt::xtensor<double, 2> elevation;
	for (int s : m_size)
	{
		dbg_in("results/mountain/bedrock-", s, elevation);

		auto holes_id = holes(elevation.shape());

		int num_holes = 0;
		for (double hp : h_prop)
		{
			int target_holes = (int)(hp *  double(holes_id.size()));

			for (int i = num_holes; i < target_holes; ++i)
				dig_hole(elevation, holes_id[i]);
			num_holes = target_holes;

			std::cerr << s << ' ' << num_holes << std::endl;
		}
	}

}