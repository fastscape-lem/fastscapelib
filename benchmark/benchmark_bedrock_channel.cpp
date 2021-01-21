#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>

#include <benchmark/benchmark.h>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"

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
        xt::xtensor<index_t, 1> receivers;
        xt::xtensor<double, 1> dist2receivers;
        xt::xtensor<index_t, 1> stack;
        xt::xtensor<double, Nd> drainage_area;

        if (Nd == 1)
        {
            double spacing = 300.;
            double x0 = 300.;
            double length = (ns - 1) * spacing;

            xt::xtensor<double, Nd> x = xt::linspace<double>(length + x0, x0, ns);
            erosion = xt::zeros_like(x);
            elevation = (length + x0 - x) * 1e-4;
            receivers = xt::arange<index_t>(-1, ns - 1);
            dist2receivers = xt::full_like(elevation, spacing);
            stack = xt::arange<index_t>(0, ns);
            drainage_area = 6.69 * xt::pow(x, 1.67);   // Hack's law
        }

        else if (Nd == 2)
        {
            auto s = bms::FastscapeSetupBase<bms::surface_type::cone, T>(state.range(0));

            fs::compute_receivers_d8(s.receivers, s.dist2receivers,
                                     s.elevation, s.active_nodes,
                                     s.dx, s.dy);
            fs::compute_donors(s.ndonors, s.donors, s.receivers);
            fs::compute_stack(s.stack, s.ndonors,
                              s.donors, s.receivers);

            erosion = xt::zeros_like(s.elevation);
            elevation = s.elevation;
            receivers = s.receivers;
            dist2receivers = s.dist2receivers;
            stack = s.stack;
            drainage_area = xt::empty_like(s.elevation);
            fs::compute_drainage_area(drainage_area, stack, receivers, s.dx, s.dy);
        }

        K k_coef = get_k_coef<K>(elevation.shape());
        double m_exp = 0.5;
        double n_exp = 1;

        double dt = 1e4;

        index_t ncorr;

        for (auto _ : state)
        {
            ncorr = fs::erode_stream_power(erosion, elevation, stack,
                                           receivers, dist2receivers, drainage_area,
                                           k_coef, m_exp, n_exp, dt, 1e-3);
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
