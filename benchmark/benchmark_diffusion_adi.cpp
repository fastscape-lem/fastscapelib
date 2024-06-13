#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>

#include <benchmark/benchmark.h>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"

#include "fastscapelib/eroders/diffusion_adi.hpp"

#include "benchmark_setup.hpp"


namespace fs = fastscapelib;
namespace bms = fs::bench_setup;


namespace fastscapelib
{

    namespace benchmark_hillslope
    {

        template <class K, class S>
        typename std::enable_if_t<std::is_floating_point<K>::value, K> get_k_coef(S&& shape)
        {
            (void) shape;  // do not show warning
            return 1e-3;
        }


        template <class K, class S>
        typename std::enable_if_t<xt::is_xexpression<K>::value, K> get_k_coef(S&& shape)
        {
            using T = typename K::value_type;

            K k_coef_arr = xt::empty<T>(shape);
            k_coef_arr.fill(1e-3);

            return k_coef_arr;
        }

        template <class T, class K>
        void bm_erode_linear_diffusion(benchmark::State& state)
        {
            auto ns = static_cast<std::size_t>(state.range(0));
            std::array<std::size_t, 2> shape = { ns, ns };

            auto grid = fs::raster_grid(shape, { 0.4, 0.4 }, fs::node_status::fixed_value);

            xt::xtensor<double, 2> elevation = xt::random::rand<T>(shape);

            auto eroder = fs::make_diffusion_adi_eroder(grid, 1e-3);
            eroder.set_k_coef(get_k_coef<K>(elevation.shape()));

            double dt = 1e4;


            for (auto _ : state)
            {
                eroder.erode(elevation, dt);
            }
        }


        BENCHMARK_TEMPLATE(bm_erode_linear_diffusion, double, double)
            ->Apply(bms::grid_sizes<benchmark::kMillisecond>);

        BENCHMARK_TEMPLATE(bm_erode_linear_diffusion, double, xt::xtensor<double, 2>)
            ->Apply(bms::grid_sizes<benchmark::kMillisecond>);

    }  // namespace benchmark_hillslope

}  // namespace fastscapelib
