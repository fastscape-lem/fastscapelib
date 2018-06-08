#include <array>

#include <benchmark/benchmark.h>

#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"

#include "fastscapelib/sinks.hpp"


namespace fastscapelib
{

namespace benchmark_sinks
{

    /*
     * Create a DEM of a nearly flat surface with small random
     * perturbation on a n x n square grid.
     *
     * The surface is characterized by a lot of small depressions.
     */
    template<class T, class S>
    auto init_random_elevation(S n) -> xt::xtensor<T, 2>
    {
        auto ns = static_cast<size_t>(n);
        std::array<size_t, 2> shape = { ns, ns };

        xt::random::seed(0);
        xt::xtensor<T, 2> elev = xt::random::rand<T>(shape);

        return elev;
    }


    /*
     * Create a DEM with "jail" walls on the boundaries of a n x n
     * square grid.
     *
     * The surface is characterized by a single, big closed depression
     * covering the whole domain (except the boundaries).
     */
    template<class T, class S>
    auto init_jail_elevation(S n) -> xt::xtensor<T, 2>
    {
        auto ns = static_cast<size_t>(n);
        std::array<size_t, 2> shape = { ns, ns };

        //TODO
    }


    static void fill_sinks_flat(benchmark::State& state)
    {
        auto elev = init_random_elevation<double>(state.range(0));

        for (auto _ : state)
        {
            fastscapelib::fill_sinks_flat(elev);
        }
    }


    static void fill_sinks_sloped(benchmark::State& state)
    {
        auto elev = init_random_elevation<double>(state.range(0));

        for (auto _ : state)
        {
            fastscapelib::fill_sinks_sloped(elev);
        }
    }


    BENCHMARK(fill_sinks_flat)
    ->Arg(100)->Arg(200)->Arg(500)->Arg(1000)->Arg(2000)->Arg(5000)
    ->Unit(benchmark::kMillisecond);

    BENCHMARK(fill_sinks_sloped)
    ->Arg(100)->Arg(200)->Arg(500)->Arg(1000)->Arg(2000)->Arg(5000)
    ->Unit(benchmark::kMillisecond);

} // namespace benchmark_sinks

} // namespace fastscapelib
