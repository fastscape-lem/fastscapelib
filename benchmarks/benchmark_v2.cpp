#include "random_benchmark.hpp"

#include "fastscapelib/flow_routing_v2.hpp"

namespace fs = fastscapelib;

void benchmark_v2(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>& active_nodes)
{
    const auto elev_shape = elevation.shape();
    const size_t nrows = (size_t) elev_shape[0];
    const size_t ncols = (size_t) elev_shape[1];
    const std::array<std::size_t, 1> shape_1d{nrows*ncols};
    xt::xtensor<size_t, 1> receivers (shape_1d);
    xt::xtensor<size_t, 1> passes (shape_1d);
    fs::compute_receivers_passes(receivers, passes, elevation, active_nodes, 1.0, 1.0);

}

//RandomBenchmark::Register register_v2("v2", benchmark_v2);
