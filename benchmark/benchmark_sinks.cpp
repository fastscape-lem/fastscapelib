#include <array>
#include <math.h>

#include <benchmark/benchmark.h>

#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"
#include "fastscapelib/sinks.hpp"

#ifdef ENABLE_RICHDEM
#include "fastscapelib/richdem.hpp"
#endif

#include "fixtures.hpp"


namespace fastscapelib
{

namespace benchmark_sinks
{


    template<fixtures::SurfaceType surf_type, class T>
    inline auto fill_sinks_flat(benchmark::State& state)
    {
        auto topo = fixtures::SyntheticTopography<surf_type, T>(state.range(0));
        auto elev = topo.get_elevation();

        for (auto _ : state)
        {
            fastscapelib::fill_sinks_flat(elev);
        }
    }


    template<fixtures::SurfaceType surf_type, class T>
    static void fill_sinks_sloped(benchmark::State& state)
    {
        auto topo = fixtures::SyntheticTopography<surf_type, T>(state.range(0));
        auto elev = topo.get_elevation();

        for (auto _ : state)
        {
            fastscapelib::fill_sinks_sloped(elev);
        }
    }

#ifdef ENABLE_RICHDEM
    template<fixtures::SurfaceType surf_type, class T>
    inline auto fill_sinks_wei2018(benchmark::State& state)
    {
        auto topo = fixtures::SyntheticTopography<surf_type, T>(state.range(0));
        auto elev = topo.get_elevation();

        for (auto _ : state)
        {
            fastscapelib::fill_sinks_wei2018(elev);
        }
    }
#endif

#ifdef ENABLE_RICHDEM
    BENCHMARK_TEMPLATE(fill_sinks_wei2018, fixtures::SurfaceType::cone, double)
    ->Arg(100)->Arg(200)->Arg(500)->Arg(1000)->Arg(2000)->Arg(5000)
    ->Unit(benchmark::kMillisecond);

    BENCHMARK_TEMPLATE(fill_sinks_wei2018, fixtures::SurfaceType::cone_inv, double)
    ->Arg(100)->Arg(200)->Arg(500)->Arg(1000)->Arg(2000)->Arg(5000)
    ->Unit(benchmark::kMillisecond);
#endif

    BENCHMARK_TEMPLATE(fill_sinks_flat, fixtures::SurfaceType::cone, double)
    ->Arg(100)->Arg(200)->Arg(500)->Arg(1000)->Arg(2000)->Arg(5000)
    ->Unit(benchmark::kMillisecond);

    BENCHMARK_TEMPLATE(fill_sinks_flat, fixtures::SurfaceType::cone_inv, double)
    ->Arg(100)->Arg(200)->Arg(500)->Arg(1000)->Arg(2000)->Arg(5000)
    ->Unit(benchmark::kMillisecond);

    BENCHMARK_TEMPLATE(fill_sinks_flat, fixtures::SurfaceType::cone_noise, double)
    ->Arg(100)->Arg(200)->Arg(500)->Arg(1000)->Arg(2000)->Arg(5000)
    ->Unit(benchmark::kMillisecond);

    BENCHMARK_TEMPLATE(fill_sinks_sloped, fixtures::SurfaceType::cone, double)
    ->Arg(100)->Arg(200)->Arg(500)->Arg(1000)->Arg(2000)->Arg(5000)
    ->Unit(benchmark::kMillisecond);

    BENCHMARK_TEMPLATE(fill_sinks_sloped, fixtures::SurfaceType::cone_inv, double)
    ->Arg(100)->Arg(200)->Arg(500)->Arg(1000)->Arg(2000)->Arg(5000)
    ->Unit(benchmark::kMillisecond);

    BENCHMARK_TEMPLATE(fill_sinks_sloped, fixtures::SurfaceType::cone_noise, double)
    ->Arg(100)->Arg(200)->Arg(500)->Arg(1000)->Arg(2000)->Arg(5000)
    ->Unit(benchmark::kMillisecond);

} // namespace benchmark_sinks

} // namespace fastscapelib
