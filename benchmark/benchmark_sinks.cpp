#include <array>
#include <math.h>

#include <benchmark/benchmark.h>

#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"

#include "fastscapelib/sinks.hpp"


namespace fastscapelib
{

namespace benchmark_sinks
{

    /*
     * A DEM (n x n square grid) representing a conic surface.
     *
     * It is a smooth surface with regular slope. It has no depression.
     */
    template<class T>
    struct conic_surface
    {
        xt::xtensor<T, 2> operator()(int n)
        {
            auto grid = xt::meshgrid(xt::linspace<double>(-1, 1, n),
                                     xt::linspace<double>(-1, 1, n));

            xt::xtensor<T, 2> elev = std::sqrt(2.) - xt::sqrt(
                xt::pow(std::get<0>(grid), 2) +
                xt::pow(std::get<1>(grid), 2));

            return elev;
        }
    };


    /*
     * A DEM (n x n square grid) representing an inverse conic
     * surface.
     *
     * This generates a single, big closed depression.
     */
    template<class T>
    struct conic_surface_inv
    {
        xt::xtensor<T, 2> operator()(int n)
        {
            xt::xtensor<T, 2> elev = -conic_surface<T>()(n);
            return elev;
        }
    };


    /*
     * A DEM (n x n square grid) representing a conic surface
     * with random perturbations.
     *
     * The magnitude of the random perturbation is chosen arbitrarily
     * so that the generated surface has many depressions of different
     * sizes.
     */
    template<class T>
    struct conic_surface_noise
    {
        xt::xtensor<T, 2> operator()(int n)
        {
            auto ns = static_cast<size_t>(n);
            std::array<size_t, 2> shape = { ns, ns };

            xt::random::seed(0);
            xt::xtensor<T, 2> elev = (conic_surface<T>()(n) +
                                      xt::random::rand<T>(shape) * 5. / n);

            return elev;
        }
    };


    /*
     * A DEM (n x n square grid) representing a nearly flat surface
     * with small random perturbations.
     *
     * This generates a lot of small depressions.
     */
    template<class T>
    struct random_elevation
    {
        xt::xtensor<T, 2> operator()(int n)
        {
            auto ns = static_cast<size_t>(n);
            std::array<size_t, 2> shape = { ns, ns };

            xt::random::seed(0);
            xt::xtensor<T, 2> elev = xt::random::rand<T>(shape);

            return elev;
        }
    };


    /*
     * A DEM (n x n square grid) representing a gaussian surface.
     *
     * It is a smooth surface that has no depression.
     */
    template<class T>
    struct gaussian_elevation
    {
        xt::xtensor<T, 2> operator()(int n)
        {
            auto grid = xt::meshgrid(xt::linspace<double>(-1, 1, n),
                                     xt::linspace<double>(-1, 1, n));

            xt::xtensor<T, 2> elev = xt::exp(
                -(xt::pow(std::get<0>(grid), 2) / 2.
                  + xt::pow(std::get<1>(grid), 2) / 2.));

            return elev;
        }
    };


    template<class Surface>
    inline auto fill_sinks_flat(benchmark::State& state)
    {
        auto elev = Surface()(state.range(0));

        for (auto _ : state)
        {
            fastscapelib::fill_sinks_flat(elev);
        }
    }


    template<class Surface>
    static void fill_sinks_sloped(benchmark::State& state)
    {
        auto elev = Surface()(state.range(0));

        for (auto _ : state)
        {
            fastscapelib::fill_sinks_sloped(elev);
        }
    }

    BENCHMARK_TEMPLATE(fill_sinks_flat, conic_surface<double>)
    ->Arg(100)->Arg(200)->Arg(500)->Arg(1000)->Arg(2000)->Arg(5000)
    ->Unit(benchmark::kMillisecond);

    BENCHMARK_TEMPLATE(fill_sinks_flat, conic_surface_inv<double>)
    ->Arg(100)->Arg(200)->Arg(500)->Arg(1000)->Arg(2000)->Arg(5000)
    ->Unit(benchmark::kMillisecond);

    BENCHMARK_TEMPLATE(fill_sinks_flat, conic_surface_noise<double>)
    ->Arg(100)->Arg(200)->Arg(500)->Arg(1000)->Arg(2000)->Arg(5000)
    ->Unit(benchmark::kMillisecond);

    BENCHMARK_TEMPLATE(fill_sinks_sloped, conic_surface<double>)
    ->Arg(100)->Arg(200)->Arg(500)->Arg(1000)->Arg(2000)->Arg(5000)
    ->Unit(benchmark::kMillisecond);

    BENCHMARK_TEMPLATE(fill_sinks_sloped, conic_surface_inv<double>)
    ->Arg(100)->Arg(200)->Arg(500)->Arg(1000)->Arg(2000)->Arg(5000)
    ->Unit(benchmark::kMillisecond);

    BENCHMARK_TEMPLATE(fill_sinks_sloped, conic_surface_noise<double>)
    ->Arg(100)->Arg(200)->Arg(500)->Arg(1000)->Arg(2000)->Arg(5000)
    ->Unit(benchmark::kMillisecond);

} // namespace benchmark_sinks

} // namespace fastscapelib
