#include <array>
#include <functional>
#include <math.h>

#include <benchmark/benchmark.h>

#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"
#include "fastscapelib/sinks.hpp"

#ifdef ENABLE_RICHDEM
#include "fastscapelib/richdem.hpp"
#endif

#include "fixtures.hpp"

#include <iostream>


using namespace fixtures;


namespace fastscapelib
{

namespace benchmark_sinks
{

    enum class Method {flat, sloped, wei2018};


    static void grid_sizes(benchmark::internal::Benchmark* bench) {
        std::array<int, 4> sizes {500, 1000, 2000, 5000};

        for (int s : sizes)
        {
            bench->Args({s})->Unit(benchmark::kMillisecond);
        }
    }


    template<Method sink_func, SurfaceType surf_type, class T>
    inline auto bm_sinks(benchmark::State& state)
    {
        auto topo = SyntheticTopography<surf_type, T>(state.range(0));
        auto elevation = topo.get_elevation();

        std::function<void(void)> fill_sinks;

        switch (sink_func) {
        case Method::flat:
            fill_sinks = [&](){ fastscapelib::fill_sinks_flat(elevation); };
            break;

        case Method::sloped:
            fill_sinks = [&](){ fastscapelib::fill_sinks_sloped(elevation); };
            break;

#ifdef ENABLE_RICHDEM
        case Method::wei2018:
            fill_sinks = [&](){ fastscapelib::fill_sinks_wei2018(elevation); };
            break;
#endif
        }

        for (auto _ : state)
        {
            fill_sinks();
        }
    }


#ifdef ENABLE_RICHDEM
    BENCHMARK_TEMPLATE(bm_sinks, Method::wei2018, SurfaceType::cone, double)
    ->Apply(grid_sizes);

    BENCHMARK_TEMPLATE(bm_sinks, Method::wei2018, SurfaceType::cone_inv, double)
    ->Apply(grid_sizes);

    BENCHMARK_TEMPLATE(bm_sinks, Method::wei2018, SurfaceType::cone_noise, double)
    ->Apply(grid_sizes);
#endif


    BENCHMARK_TEMPLATE(bm_sinks, Method::flat, SurfaceType::cone, double)
    ->Apply(grid_sizes);

    BENCHMARK_TEMPLATE(bm_sinks, Method::flat, SurfaceType::cone_inv, double)
    ->Apply(grid_sizes);

    BENCHMARK_TEMPLATE(bm_sinks, Method::flat, SurfaceType::cone_noise, double)
    ->Apply(grid_sizes);

    BENCHMARK_TEMPLATE(bm_sinks, Method::sloped, SurfaceType::cone, double)
    ->Apply(grid_sizes);

    BENCHMARK_TEMPLATE(bm_sinks, Method::sloped, SurfaceType::cone_inv, double)
    ->Apply(grid_sizes);

    BENCHMARK_TEMPLATE(bm_sinks, Method::sloped, SurfaceType::cone_noise, double)
    ->Apply(grid_sizes);

} // namespace benchmark_sinks

} // namespace fastscapelib
