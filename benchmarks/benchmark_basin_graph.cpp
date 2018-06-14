#include "random_benchmark.hpp"
#include "fastscape_benchmark.hpp"

#include "fastscapelib/flow_routing.hpp"

void benchmark_basin_flat_kruskal_simple(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>& active_nodes)
{

    fs::fill_sinks_flat_basin_graph<fs::BasinAlgo::Kruskal, fs::ConnectType::Simple>(elevation, active_nodes, 1.0, 1.0);

}

RegisterRandom register_basin_flat_kruskal_simple("Kruskal simple", benchmark_basin_flat_kruskal_simple);

void benchmark_basin_flat_boruvka_simple(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>& active_nodes)
{

    fs::fill_sinks_flat_basin_graph<fs::BasinAlgo::Boruvka, fs::ConnectType::Simple>(elevation, active_nodes, 1.0, 1.0);

}

RegisterRandom register_basin_flat_boruvka("Boruvka simple", benchmark_basin_flat_boruvka_simple);

void benchmark_basin_flat_kruskal_carved(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>& active_nodes)
{

    fs::fill_sinks_flat_basin_graph<fs::BasinAlgo::Kruskal, fs::ConnectType::Carved>(elevation, active_nodes, 1.0, 1.0);

}

RegisterRandom register_basin_flat_kruskal_carved("Kruskal carved", benchmark_basin_flat_kruskal_carved);

void benchmark_basin_flat_boruvka_carved(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>& active_nodes)
{

    fs::fill_sinks_flat_basin_graph<fs::BasinAlgo::Boruvka, fs::ConnectType::Carved>(elevation, active_nodes, 1.0, 1.0);

}

RegisterRandom register_basin_flat_boruvka_carved("Boruvka carved", benchmark_basin_flat_boruvka_carved);

void benchmark_basin_flat_kruskal_sloped(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>& active_nodes)
{

    fs::fill_sinks_flat_basin_graph<fs::BasinAlgo::Kruskal, fs::ConnectType::Sloped>(elevation, active_nodes, 1.0, 1.0);

}

RegisterRandom register_basin_flat_kruskal_sloped("Kruskal sloped", benchmark_basin_flat_kruskal_sloped);

void benchmark_basin_flat_boruvka_sloped(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>& active_nodes)
{

    fs::fill_sinks_flat_basin_graph<fs::BasinAlgo::Boruvka, fs::ConnectType::Sloped>(elevation, active_nodes, 1.0, 1.0);

}

RegisterRandom register_basin_flat_boruvka_sloped("Boruvka sloped", benchmark_basin_flat_boruvka_sloped);

template<fs::BasinAlgo algo, fs::ConnectType connect>
void benchmark_fastscape_basin(
	xt::xtensor<index_t, 1>&      stack,
	xt::xtensor<index_t, 1>&      receivers,
	xt::xtensor<double, 1>&       dist2receivers,
	const xt::xtensor<double, 2>& elevation,
	const xt::xtensor<bool, 2>&   active_nodes,
	double dx, double dy)
{
	const auto elev_shape = elevation.shape();
	const size_t nrows = (size_t)elev_shape[0];
	const size_t ncols = (size_t)elev_shape[1];
	std::array<std::size_t, 1> shape_1D{ nrows * ncols };

	fs::BasinGraph<index_t, index_t, double> basin_graph;

	xt::xtensor<index_t, 1> basins(shape_1D);
	xt::xtensor<index_t, 1> ndonors(shape_1D);
	xt::xtensor<index_t, 2> donors({ nrows * ncols, 8 });


	fs::compute_receivers_d8(receivers, dist2receivers,
		elevation, active_nodes,
		dx, dy);

	fs::compute_donors(ndonors, donors, receivers);
	fs::compute_stack(stack, ndonors, donors, receivers);

	fs::correct_flowrouting<algo, connect>(basin_graph, basins, receivers, dist2receivers,
		ndonors, donors, stack, active_nodes, elevation, dx, dy);
}

//RegisterFastscape register_fastscape_basin_ks("Kruskal simple", benchmark_fastscape_basin<fs::BasinAlgo::Kruskal, fs::ConnectType::Simple>);
//RegisterFastscape register_fastscape_basin_bs("Boruvka simple", benchmark_fastscape_basin<fs::BasinAlgo::Boruvka, fs::ConnectType::Simple>);
//RegisterFastscape register_fastscape_basin_kc("Kruskal carved", benchmark_fastscape_basin<fs::BasinAlgo::Kruskal, fs::ConnectType::Carved>);
//RegisterFastscape register_fastscape_basin_bc("Boruvka carved", benchmark_fastscape_basin<fs::BasinAlgo::Boruvka, fs::ConnectType::Carved>);
RegisterFastscape register_fastscape_basin_ksl("Kruskal sloped", benchmark_fastscape_basin<fs::BasinAlgo::Kruskal, fs::ConnectType::Sloped>);
//RegisterFastscape register_fastscape_basin_bsl("Boruvka sloped", benchmark_fastscape_basin<fs::BasinAlgo::Boruvka, fs::ConnectType::Sloped>);
