#include "random_benchmark.hpp"
#include "fastscape_benchmark.hpp"

#include "fastscapelib/flow_routing.hpp"

#include "benchmark_basin_graph.hpp"


namespace fs = fastscapelib;

/*
void benchmark_basin_flat_kruskal_simple(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>&
active_nodes)
{

    fs::fill_sinks_flat_basin_graph<fs::BasinAlgo::Kruskal, fs::ConnectType::Simple>(elevation,
active_nodes, 1.0, 1.0);

}

RegisterRandom register_basin_flat_kruskal_simple("Kruskal simple",
benchmark_basin_flat_kruskal_simple);

void benchmark_basin_flat_boruvka_simple(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>&
active_nodes)
{

    fs::fill_sinks_flat_basin_graph<fs::BasinAlgo::Boruvka, fs::ConnectType::Simple>(elevation,
active_nodes, 1.0, 1.0);

}

RegisterRandom register_basin_flat_boruvka("Boruvka simple", benchmark_basin_flat_boruvka_simple);

void benchmark_basin_flat_kruskal_carved(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>&
active_nodes)
{

    fs::fill_sinks_flat_basin_graph<fs::BasinAlgo::Kruskal, fs::ConnectType::Carved>(elevation,
active_nodes, 1.0, 1.0);

}

RegisterRandom register_basin_flat_kruskal_carved("Kruskal carved",
benchmark_basin_flat_kruskal_carved);

void benchmark_basin_flat_boruvka_carved(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>&
active_nodes)
{

    fs::fill_sinks_flat_basin_graph<fs::BasinAlgo::Boruvka, fs::ConnectType::Carved>(elevation,
active_nodes, 1.0, 1.0);

}

RegisterRandom register_basin_flat_boruvka_carved("Boruvka carved",
benchmark_basin_flat_boruvka_carved);

void benchmark_basin_flat_kruskal_sloped(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>&
active_nodes)
{

    fs::fill_sinks_flat_basin_graph<fs::BasinAlgo::Kruskal, fs::ConnectType::Sloped>(elevation,
active_nodes, 1.0, 1.0);

}

RegisterRandom register_basin_flat_kruskal_sloped("Kruskal sloped",
benchmark_basin_flat_kruskal_sloped);

void benchmark_basin_flat_boruvka_sloped(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>&
active_nodes)
{

    fs::fill_sinks_flat_basin_graph<fs::BasinAlgo::Boruvka, fs::ConnectType::Sloped>(elevation,
active_nodes, 1.0, 1.0);

}

RegisterRandom register_basin_flat_boruvka_sloped("Boruvka sloped",
benchmark_basin_flat_boruvka_sloped);


//RegisterFastscape register_fastscape_basin_ks("Kruskal simple",
benchmark_fastscape_basin<fs::BasinAlgo::Kruskal, fs::ConnectType::Simple>);
//RegisterFastscape register_fastscape_basin_bs("Boruvka simple",
benchmark_fastscape_basin<fs::BasinAlgo::Boruvka, fs::ConnectType::Simple>); RegisterFastscape
register_fastscape_basin_kc("Kruskal carved", benchmark_fastscape_basin<fs::BasinAlgo::Kruskal,
fs::ConnectType::Carved>); RegisterFastscape register_fastscape_basin_bc("Boruvka carved",
benchmark_fastscape_basin<fs::BasinAlgo::Boruvka, fs::ConnectType::Carved>); RegisterFastscape
register_fastscape_basin_ksl("Kruskal sloped", benchmark_fastscape_basin<fs::BasinAlgo::Kruskal,
fs::ConnectType::Sloped>); RegisterFastscape register_fastscape_basin_bsl("Boruvka sloped",
benchmark_fastscape_basin<fs::BasinAlgo::Boruvka, fs::ConnectType::Sloped>);
*/
