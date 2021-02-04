#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>

#include <benchmark/benchmark.h>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"

#include "fastscapelib/raster_grid.hpp"
#include "fastscapelib/flow_graph.hpp"
#include "fastscapelib/flow_routing.hpp"
#include "fastscapelib/bedrock_channel.hpp"

#include "benchmark_setup.hpp"


namespace fs = fastscapelib;
namespace bms = fs::bench_setup;


namespace fastscapelib
{
    namespace benchmark_bedrock_channel
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


        template<class T, class K, int Nd>
        void bm_erode_stream_power(benchmark::State& state)
        {
            auto ns = static_cast<std::size_t>(state.range(0));

            xt::xtensor<double, Nd> erosion;
            xt::xtensor<double, Nd> elevation;

            if (Nd == 1)
            {
                using flow_graph_type = fs::flow_graph<fs::profile_grid, double>;
                using index_type = typename flow_graph_type::index_type;

                double spacing = 300.;

                auto grid = fs::profile_grid(ns, spacing, fs::node_status::fixed_value_boundary);
            
                auto flow_graph = flow_graph_type(grid,
                                                  std::make_unique<fs::single_flow_router<flow_graph_type>>(),
                                                  std::make_unique<fs::no_sink_resolver<flow_graph_type>>());

                double x0 = 300.;
                double length = (ns - 1) * spacing;

                xt::xtensor<double, 1> x = xt::linspace<double>(length + x0, x0, ns);
                erosion = xt::zeros_like(x);
                elevation = (length + x0 - x) * 1e-4;

                flow_graph.update_routes(elevation);

                auto drainage_area = flow_graph.accumulate(1.);

                K k_coef = get_k_coef<K>(elevation.shape());
                double m_exp = 0.5;
                double n_exp = 1;
                double dt = 1e4;

                index_type ncorr;

                for (auto _ : state)
                {
                    ncorr = fs::erode_stream_power(erosion, elevation, drainage_area, flow_graph,
                                                   k_coef, m_exp, n_exp, dt, 1e-3);
                }
            }

            else if (Nd == 2)
            {
                using flow_graph_type = fs::flow_graph<fs::raster_grid, double>;
                using index_type = typename flow_graph_type::index_type;

                auto s = bms::FastscapeSetupBase<bms::surface_type::cone, T>(state.range(0));

                auto grid = fs::raster_grid({{ns, ns}}, {s.dy, s.dy}, fs::node_status::fixed_value_boundary);
            
                auto flow_graph = flow_graph_type(grid,
                                                  std::make_unique<fs::single_flow_router<flow_graph_type>>(),
                                                  std::make_unique<fs::no_sink_resolver<flow_graph_type>>());

                elevation = s.elevation;
                erosion = xt::zeros_like(s.elevation);
                flow_graph.update_routes(elevation);

                auto drainage_area = flow_graph.accumulate(1.);

                K k_coef = get_k_coef<K>(elevation.shape());
                double m_exp = 0.5;
                double n_exp = 1;
                double dt = 1e4;

                index_type ncorr;

                for (auto _ : state)
                {
                    ncorr = fs::erode_stream_power(erosion, elevation, drainage_area, flow_graph,
                                                   k_coef, m_exp, n_exp, dt, 1e-3);
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
