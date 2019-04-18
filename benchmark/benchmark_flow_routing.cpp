#include <benchmark/benchmark.h>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xtensor.hpp"

#include "fastscapelib/flow_routing.hpp"

#include "benchmark_setup.hpp"


namespace fs = fastscapelib;
namespace bms = benchmark_setup;


namespace fastscapelib
{

namespace benchmark_grid
{

template<class T>
void bm_compute_receivers_d8(benchmark::State& state)
{
    auto s = bms::FastscapeSetupBase<bms::SurfaceType::flat_noise, T>(state.range(0));

    for (auto _ : state)
    {
        fs::compute_receivers_d8(s.receivers, s.dist2receivers,
                                 s.elevation, s.active_nodes,
                                 s.dx, s.dy);
    }
}

BENCHMARK_TEMPLATE(bm_compute_receivers_d8, double)
->Apply(bms::grid_sizes<benchmark::kMillisecond>);

} // namespace benchmark_grid

} // namespace fastscapelib
