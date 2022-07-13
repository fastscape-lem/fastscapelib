#include <benchmark/benchmark.h>

#include "fastscapelib/algo/pflood.hpp"

#include "benchmark_setup.hpp"


namespace fs = fastscapelib;
namespace bms = fs::bench_setup;


namespace fastscapelib
{
    namespace bench
    {

        using surf_type = bms::surface_type;


        template <surf_type S>
        void fill_sinks_flat(benchmark::State& state)
        {
            auto topo = bms::synthetic_topography_2d<S, double>(state.range(0));
            auto grid = topo.grid();
            auto elevation = topo.elevation();

            for (auto _ : state)
            {
                fs::fill_sinks_flat(grid, elevation);
            }
        }

        template <surf_type S>
        void fill_sinks_sloped(benchmark::State& state)
        {
            auto topo = bms::synthetic_topography_2d<S, double>(state.range(0));
            auto grid = topo.grid();
            auto elevation = topo.elevation();

            for (auto _ : state)
            {
                fs::fill_sinks_sloped(grid, elevation);
            }
        }


#define BENCH_SURF(NAME, SURF)                                                                     \
    BENCHMARK_TEMPLATE(NAME, SURF)->Apply(bms::small_grid_sizes<benchmark::kMillisecond>);

#define BENCH_ALL_SURF(NAME)                                                                       \
    BENCH_SURF(NAME, surf_type::cone)                                                              \
    BENCH_SURF(NAME, surf_type::cone_inv)                                                          \
    BENCH_SURF(NAME, surf_type::cone_noise)                                                        \
    BENCH_SURF(NAME, surf_type::flat_noise)


        BENCH_ALL_SURF(fill_sinks_flat);
        BENCH_ALL_SURF(fill_sinks_sloped);

    }
}
