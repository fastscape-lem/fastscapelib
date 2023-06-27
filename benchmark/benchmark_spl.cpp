#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>

#include <benchmark/benchmark.h>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"

#include "fastscapelib/grid/raster_grid.hpp"
#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/eroders/spl.hpp"

#include "benchmark_setup.hpp"


namespace fs = fastscapelib;
namespace bms = fs::bench_setup;


namespace fastscapelib
{
    namespace benchmark_bedrock_channel
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


        template <class T, class K, int Nd>
        void bm_erode_stream_power(benchmark::State& state)
        {
            auto ns = static_cast<std::size_t>(state.range(0));

            xt::xtensor<double, Nd> erosion;
            xt::xtensor<double, Nd> elevation;

            if (Nd == 1)
            {
                using flow_graph_type = fs::flow_graph<fs::profile_grid>;

                double spacing = 300.;
                double x0 = 300.;
                double length = (ns - 1.0) * spacing;
                xt::xtensor<double, 1> x = xt::linspace<double>(length + x0, x0, ns);
                elevation = (length + x0 - x) * 1e-4;

                auto grid = fs::profile_grid(ns, spacing, fs::node_status::fixed_value);

                auto flow_graph = flow_graph_type(grid, { fs::single_flow_router() });

                K k_coef = get_k_coef<K>(elevation.shape());
                double area_exp = 0.5;
                double slope_exp = 1;
                auto eroder = fs::make_spl_eroder(flow_graph, k_coef, area_exp, slope_exp, 1e-3);

                flow_graph.update_routes(elevation);
                auto drainage_area = flow_graph.accumulate(1.);

                double dt = 1e4;

                for (auto _ : state)
                {
                    eroder.erode(elevation, drainage_area, dt);
                }
            }

            else if (Nd == 2)
            {
                using flow_graph_type = fs::flow_graph<fs::raster_grid>;

                auto s = bms::FastscapeSetupBase<bms::surface_type::cone, T>(state.range(0));
                elevation = s.elevation;

                auto grid
                    = fs::raster_grid({ { ns, ns } }, { s.dy, s.dy }, fs::node_status::fixed_value);

                auto flow_graph = flow_graph_type(grid, { fs::single_flow_router() });

                K k_coef = get_k_coef<K>(elevation.shape());
                double area_exp = 0.5;
                double slope_exp = 1;
                auto eroder = fs::make_spl_eroder(flow_graph, k_coef, area_exp, slope_exp, 1e-3);

                flow_graph.update_routes(elevation);
                auto drainage_area = flow_graph.accumulate(1.);

                double dt = 1e4;

                for (auto _ : state)
                {
                    eroder.erode(elevation, drainage_area, dt);
                }
            }
        }


        BENCHMARK_TEMPLATE(bm_erode_stream_power, double, double, 1)
            ->Apply(bms::grid_sizes<benchmark::kMicrosecond>);

        BENCHMARK_TEMPLATE(bm_erode_stream_power, double, xt::xtensor<double, 1>, 1)
            ->Apply(bms::grid_sizes<benchmark::kMicrosecond>);

        BENCHMARK_TEMPLATE(bm_erode_stream_power, double, double, 2)
            ->Apply(bms::grid_sizes<benchmark::kMillisecond>);

        BENCHMARK_TEMPLATE(bm_erode_stream_power, double, xt::xtensor<double, 2>, 2)
            ->Apply(bms::grid_sizes<benchmark::kMillisecond>);

    }  // namespace benchmark_bedrock_channel
}  // namespace fastscapelib
