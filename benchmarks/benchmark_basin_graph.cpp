#include "random_benchmark.hpp"

#include "fastscapelib/flow_routing.hpp"

void benchmark_basin_flat_kruskal(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>& active_nodes)
{

    fs::fill_sinks_flat_basin_graph<fs::BasinAlgo::Kruskal>(elevation, active_nodes, 1.0, 1.0);

}

RandomBenchmark::Register register_basin_flat_kruskal("Basin Flat kruskal", benchmark_basin_flat_kruskal);

void benchmark_basin_flat_boruvka(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>& active_nodes)
{

    fs::fill_sinks_flat_basin_graph<fs::BasinAlgo::Boruvka>(elevation, active_nodes, 1.0, 1.0);

}

RandomBenchmark::Register register_basin_flat_boruvka("Basin Flat boruvka", benchmark_basin_flat_boruvka);
