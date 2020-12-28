#include "benchmark_setup.hpp"

#include "fastscapelib/consts.hpp"
#include "fastscapelib/profile_grid.hpp"
#include "fastscapelib/raster_grid.hpp"
#include "fastscapelib/flow_routing.hpp"

#include "xtensor/xfixed.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"

#include <benchmark/benchmark.h>


namespace fs = fastscapelib;
namespace bms = benchmark_setup;


namespace fastscapelib
{
    namespace bench
    {

        template <class G>
        void flow_routing__compute_receivers__profile(benchmark::State& state)
        {
            using grid = G;
            using size_type = typename grid::size_type;

            auto n = static_cast<size_type>(state.range(0));
            auto profile_grid = grid(n, 1., fs::node_status::fixed_value_boundary);

            xt::xtensor<double, 1> elevation = xt::random::rand<double>({n});
            xt::xtensor<index_t, 1> receivers = xt::ones<index_t>({n}) * -1;
            xt::xtensor<double, 1> dist2receivers = xt::ones<double>({n}) * -1.;

            //warm-up cache
            profile_grid.neighbors_cache().warm_up();

            for (auto _ : state)
            {
                fs::compute_receivers(receivers, dist2receivers,
                                      elevation, profile_grid); 
            }
        }

        template <class G>
        void flow_routing__compute_receivers__raster(benchmark::State& state)
        {
            using grid = G;
            using size_type = typename grid::size_type;

            auto n = static_cast<size_type>(state.range(0));
            std::array<size_type, 2> shape {n, n};
            auto raster_grid = grid(shape, {1., 1.}, fs::node_status::fixed_value_boundary);

            xt::xtensor<double, 2> elevation = xt::random::rand<double>({n,n});
            xt::xtensor<index_t, 1> receivers = xt::ones<index_t>({n*n}) * -1;
            xt::xtensor<double, 1> dist2receivers = xt::ones<double>({n*n}) * -1.;

            //warm-up cache
            for (size_type idx = 0; idx < raster_grid.size(); ++idx)
            {
                raster_grid.neighbors(idx);
            }  

            for (auto _ : state)
            {
                fs::compute_receivers(receivers, dist2receivers,
                                      elevation, raster_grid); 
            }
        }

        void flow_routing__compute_receivers_d8(benchmark::State& state)
        {
            auto n = static_cast<std::size_t>(state.range(0));

            xt::xtensor<double, 2> elevation = xt::random::rand<double>({n,n});
            xt::xtensor<index_t, 1> receivers = xt::ones<index_t>({n*n}) * -1;
            xt::xtensor<double, 1> dist2receivers = xt::ones<double>({n*n}) * -1.;
            xt::xtensor<bool, 2> active_nodes = xt::full_like(elevation, true);           

            for (auto _ : state)
            {
                fs::compute_receivers_d8(receivers, dist2receivers,
                                         elevation, active_nodes, 1., 1.); 
            }
        }

        BENCHMARK_TEMPLATE(flow_routing__compute_receivers__profile, fs::profile_grid)
        ->Apply(bms::grid_sizes<benchmark::kMicrosecond>);

        BENCHMARK_TEMPLATE(flow_routing__compute_receivers__raster, fs::raster_grid)
        ->Apply(bms::grid_sizes<benchmark::kMillisecond>)->MinTime(2.);

        BENCHMARK(flow_routing__compute_receivers_d8)
        ->Apply(bms::grid_sizes<benchmark::kMillisecond>)->MinTime(2.);

    } // namespace bench
} // namespace fastscapelib