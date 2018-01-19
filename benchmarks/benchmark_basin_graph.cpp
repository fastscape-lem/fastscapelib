#include "random_benchmark.hpp"

#include "fastscapelib/flow_routing.hpp"

void benchmark_basin_flat(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>& active_nodes)
{

    fs::fill_sinks_flat_basin_graph(elevation, active_nodes, 1.0, 1.0);

}

RandomBenchmark::Register register_basin_flat("Basin Flat", benchmark_basin_flat);
