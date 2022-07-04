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

#include "benchmark_setup.hpp"


namespace bms = benchmark_setup;
namespace fs = fastscapelib;

using SurfaceType = bms::SurfaceType;


namespace fastscapelib
{

    namespace benchmark_sinks
    {

        enum class Method
        {
            flat,
            sloped,
            wei2018
        };


        template <Method sink_func, SurfaceType surf_type, class T>
        inline auto bm_sinks(benchmark::State& state)
        {
            auto topo = bms::SyntheticTopography<surf_type, T>(state.range(0));
            auto elevation = topo.get_elevation();

            std::function<void(void)> fill_sinks;

            switch (sink_func)
            {
                case Method::flat:
                    fill_sinks = [&]() { fs::fill_sinks_flat(elevation); };
                    break;

                case Method::sloped:
                    fill_sinks = [&]() { fs::fill_sinks_sloped(elevation); };
                    break;

#ifdef ENABLE_RICHDEM
                case Method::wei2018:
                    fill_sinks = [&]() { fs::fill_sinks_wei2018(elevation); };
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
            ->Apply(bms::grid_sizes);

        BENCHMARK_TEMPLATE(bm_sinks, Method::wei2018, SurfaceType::cone_inv, double)
            ->Apply(bms::grid_sizes);

        BENCHMARK_TEMPLATE(bm_sinks, Method::wei2018, SurfaceType::cone_noise, double)
            ->Apply(bms::grid_sizes);
#endif

        BENCHMARK_TEMPLATE(bm_sinks, Method::flat, SurfaceType::cone, double)
            ->Apply(bms::grid_sizes);

        BENCHMARK_TEMPLATE(bm_sinks, Method::flat, SurfaceType::cone_inv, double)
            ->Apply(bms::grid_sizes);

        BENCHMARK_TEMPLATE(bm_sinks, Method::flat, SurfaceType::cone_noise, double)
            ->Apply(bms::grid_sizes);

        BENCHMARK_TEMPLATE(bm_sinks, Method::sloped, SurfaceType::cone, double)
            ->Apply(bms::grid_sizes);

        BENCHMARK_TEMPLATE(bm_sinks, Method::sloped, SurfaceType::cone_inv, double)
            ->Apply(bms::grid_sizes);

        BENCHMARK_TEMPLATE(bm_sinks, Method::sloped, SurfaceType::cone_noise, double)
            ->Apply(bms::grid_sizes);

    }  // namespace benchmark_sinks

}  // namespace fastscapelib
