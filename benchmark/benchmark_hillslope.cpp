#include <array>
#include <cmath>
#include <type_traits>

#include <benchmark/benchmark.h>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"

#include "fastscapelib/hillslope.hpp"

#include "benchmark_setup.hpp"


namespace fs = fastscapelib;
namespace bms = benchmark_setup;


namespace fastscapelib
{

namespace benchmark_hillslope
{

    template<class K, class S>
    typename std::enable_if_t<std::is_floating_point<K>::value, K>
    get_k_coef(S&& shape)
    {
        (void) shape;  // do not show warning
        return 1e-3;
    }


    template<class K, class S>
    typename std::enable_if_t<xt::is_xexpression<K>::value, K>
    get_k_coef(S&& shape)
    {
        using T = typename K::value_type;

        K k_coef_arr = xt::empty<T>(shape);
        k_coef_arr.fill(1e-3);

        return k_coef_arr;
    }

    template<class T, class K>
    void bm_erode_linear_diffusion(benchmark::State& state)
    {
        auto ns = static_cast<size_t>(state.range(0));
        std::array<size_t, 2> shape = { ns, ns };

        xt::xtensor<double, 2> elevation = xt::random::rand<T>(shape);
        auto erosion = xt::empty_like(elevation);

        double dx = 0.4;
        double dy = 0.4;
        double dt = 1e4;

        K k_coef = get_k_coef<K>(elevation.shape());

        for (auto _ : state)
        {
            fs::erode_linear_diffusion(erosion, elevation,
                                       k_coef, dt, dx, dy);
        }
    }


    BENCHMARK_TEMPLATE(bm_erode_linear_diffusion, double, double)
    ->Apply(bms::grid_sizes<benchmark::kMillisecond>);

    BENCHMARK_TEMPLATE(bm_erode_linear_diffusion, double, xt::xtensor<double, 2>)
    ->Apply(bms::grid_sizes<benchmark::kMillisecond>);

}  // namespace benchmark_hillslope

}  // namespace fastscapelib
