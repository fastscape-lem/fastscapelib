#include "benchmark_setup.hpp"

#include "fastscapelib/grid/raster_grid.hpp"
#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"

#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"

#include <benchmark/benchmark.h>


namespace fs = fastscapelib;
namespace bms = fs::bench_setup;


namespace fastscapelib
{
    namespace bench
    {

        void single_flow_router__raster(benchmark::State& state)
        {
            using grid_type = fs::raster_grid;
            using size_type = typename grid_type::size_type;

            auto n = static_cast<size_type>(state.range(0));
            std::array<size_type, 2> shape{ { n, n } };
            auto grid = grid_type(shape, { 1., 1. }, fs::node_status::fixed_value);
            auto graph = fs::flow_graph<grid_type>(grid, { fs::single_flow_router() });

            xt::xtensor<double, 2> elevation = xt::random::rand<double>({ n, n });

            // warm-up grid cache
            for (size_type idx = 0; idx < grid.size(); ++idx)
            {
                grid.neighbors(idx);
            }

            for (auto _ : state)
            {
                graph.update_routes(elevation);
            }
        }

        BENCHMARK(single_flow_router__raster)->Apply(bms::grid_sizes<benchmark::kMillisecond>);

    }  // namespace bench
}  // namespace fastscapelib
