#include "benchmark_setup.hpp"

#include "fastscapelib/utils/consts.hpp"
#include "fastscapelib/grid/profile_grid.hpp"
#include "fastscapelib/grid/raster_grid.hpp"
#include "fastscapelib/algo/flow_routing.hpp"

#include "xtensor/xfixed.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"

#include <benchmark/benchmark.h>


namespace fs = fastscapelib;
namespace bms = fs::bench_setup;


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
            auto profile_grid = grid(n, 1., fs::node_status::fixed_value);

            xt::xtensor<double, 1> elevation = xt::random::rand<double>({ n });
            xt::xtensor<index_t, 1> receivers = xt::ones<index_t>({ n }) * -1;
            xt::xtensor<double, 1> dist2receivers = xt::ones<double>({ n }) * -1.;

            // warm-up cache
            for (size_type idx = 0; idx < profile_grid.size(); ++idx)
            {
                profile_grid.neighbors(idx);
            }

            for (auto _ : state)
            {
                fs::compute_receivers(receivers, dist2receivers, elevation, profile_grid);
            }
        }

        template <class G>
        void flow_routing__compute_receivers__raster(benchmark::State& state)
        {
            using grid = G;
            using size_type = typename grid::size_type;

            auto n = static_cast<size_type>(state.range(0));
            std::array<size_type, 2> shape{ { n, n } };
            auto raster_grid = grid(shape, { 1., 1. }, fs::node_status::fixed_value);

            xt::xtensor<double, 2> elevation = xt::random::rand<double>({ n, n });
            xt::xtensor<index_t, 1> receivers = xt::ones<index_t>({ n * n }) * -1;
            xt::xtensor<double, 1> dist2receivers = xt::ones<double>({ n * n }) * -1.;

            // warm-up cache
            for (size_type idx = 0; idx < raster_grid.size(); ++idx)
            {
                raster_grid.neighbors(idx);
            }

            for (auto _ : state)
            {
                fs::compute_receivers(receivers, dist2receivers, elevation, raster_grid);
            }
        }

        void flow_routing__compute_receivers_d8(benchmark::State& state)
        {
            auto n = static_cast<std::size_t>(state.range(0));

            xt::xtensor<double, 2> elevation = xt::random::rand<double>({ n, n });
            xt::xtensor<index_t, 1> receivers = xt::ones<index_t>({ n * n }) * -1;
            xt::xtensor<double, 1> dist2receivers = xt::ones<double>({ n * n }) * -1.;
            xt::xtensor<bool, 2> active_nodes = xt::full_like(elevation, true);

            for (auto _ : state)
            {
                fs::compute_receivers_d8(
                    receivers, dist2receivers, elevation, active_nodes, 1., 1.);
            }
        }

        using profile_nocache = fs::profile_grid_xt<xt_selector, neighbors_no_cache<2>>;
        using profile_cacheall = fs::profile_grid;

        using queen_nocache
            = fs::raster_grid_xt<xt_selector, raster_connect::queen, neighbors_no_cache<8>>;
        using queen_cacheall
            = fs::raster_grid_xt<xt_selector, raster_connect::queen, neighbors_cache<8>>;

        using rook_nocache
            = fs::raster_grid_xt<xt_selector, raster_connect::rook, neighbors_no_cache<4>>;
        using rook_cacheall
            = fs::raster_grid_xt<xt_selector, raster_connect::rook, neighbors_cache<4>>;

        using bishop_nocache
            = fs::raster_grid_xt<xt_selector, raster_connect::bishop, neighbors_no_cache<4>>;
        using bishop_cacheall
            = fs::raster_grid_xt<xt_selector, raster_connect::bishop, neighbors_cache<4>>;


#define BENCH_GRID(NAME, GRID)                                                                     \
    BENCHMARK_TEMPLATE(NAME, GRID)->Apply(bms::grid_sizes<benchmark::kMillisecond>);

#define BENCH_RASTER(NAME)                                                                         \
    BENCH_GRID(NAME, queen_nocache)                                                                \
    BENCH_GRID(NAME, queen_cacheall)                                                               \
    BENCH_GRID(NAME, rook_nocache)                                                                 \
    BENCH_GRID(NAME, rook_cacheall)                                                                \
    BENCH_GRID(NAME, bishop_nocache)                                                               \
    BENCH_GRID(NAME, bishop_cacheall)


        BENCHMARK_TEMPLATE(flow_routing__compute_receivers__profile, profile_nocache)
            ->Apply(bms::grid_sizes<benchmark::kMicrosecond>);
        BENCHMARK_TEMPLATE(flow_routing__compute_receivers__profile, profile_cacheall)
            ->Apply(bms::grid_sizes<benchmark::kMicrosecond>);

        BENCH_RASTER(flow_routing__compute_receivers__raster)

        BENCHMARK(flow_routing__compute_receivers_d8)
            ->Apply(bms::grid_sizes<benchmark::kMillisecond>);

    }  // namespace bench
}  // namespace fastscapelib
