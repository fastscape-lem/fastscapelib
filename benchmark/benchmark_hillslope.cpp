#include <array>
#include <math.h>

#include <benchmark/benchmark.h>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"

#include "fastscapelib/hillslope.hpp"


namespace fs = fastscapelib;


namespace fastscapelib
{

namespace benchmark_hillslope
{
    enum class KCoefType {scalar, array};


    template<class T, KCoefType k_coef_type>
    void bm_erode_linear_diffusion(benchmark::State& state)
    {
        auto ns = static_cast<size_t>(state.range(0));
        std::array<size_t, 2> shape = { ns, ns };

        xt::xtensor<double, 2> elevation = xt::random::rand<T>(shape);
        auto erosion = xt::empty_like(elevation);

        double dx = 0.4;
        double dy = 0.4;
        double k_coef = 1e-3;
        double dt = 1e4;

        if (k_coef_type == KCoefType::scalar)
        {
            double k_coef = 1e-3;
        }
        else if (k_coef_type == KCoefType::array)
        {
            auto k_coef = xt::full_like(elevation, 1e-3);
        }

        for (auto _ : state)
        {
            fs::erode_linear_diffusion(erosion, elevation,
                                       k_coef, dt, dx, dy);
        }
    }


    BENCHMARK_TEMPLATE(bm_erode_linear_diffusion, double, KCoefType::scalar)
    ->Arg(256)->Arg(512)->Arg(1024)->Arg(2048)->Arg(4096)
    ->Unit(benchmark::kMillisecond);

    BENCHMARK_TEMPLATE(bm_erode_linear_diffusion, double, KCoefType::array)
    ->Arg(256)->Arg(512)->Arg(1024)->Arg(2048)->Arg(4096)
    ->Unit(benchmark::kMillisecond);

}  // namespace benchmark_hillslope

}  // namespace fastscapelib
