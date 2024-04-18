#include "benchmark_setup.hpp"

#include "fastscapelib/utils/consts.hpp"
#include "fastscapelib/grid/profile_grid.hpp"
#include "fastscapelib/grid/raster_grid.hpp"
#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"

#include "xtensor/xfixed.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"

#include <benchmark/benchmark.h>


namespace fs = fastscapelib;
namespace bms = fs::bench_setup;


namespace fastscapelib::bench
{
    template <class G>
    void flow_routing__single_flow_router__profile(benchmark::State& state)
    {
        using grid = G;
        using size_type = typename grid::size_type;

        auto n = static_cast<size_type>(state.range(0));
        auto profile_grid = grid(n, 1., fs::node_status::fixed_value);

        xt::xtensor<double, 1> elevation = xt::random::rand<double>({ n });
        flow_graph<G> graph(profile_grid, { single_flow_router() });

        // warm-up cache
        for (size_type idx = 0; idx < profile_grid.size(); ++idx)
        {
            profile_grid.neighbors(idx);
        }

        for (auto _ : state)
        {
            graph.update_routes(elevation);
        }
    }

    template <class G>
    void flow_routing__single_flow_router__raster(benchmark::State& state)
    {
        using grid = G;
        using size_type = typename grid::size_type;

        auto n = static_cast<size_type>(state.range(0));
        std::array<size_type, 2> shape{ { n, n } };
        auto raster_grid = grid(shape, { 1., 1. }, fs::node_status::fixed_value);

        xt::xtensor<double, 2> elevation = xt::random::rand<double>(raster_grid.shape());
        flow_graph<grid> graph(raster_grid, { single_flow_router() });

        // warm-up cache
        for (size_type idx = 0; idx < raster_grid.size(); ++idx)
        {
            raster_grid.neighbors(idx);
        }

        for (auto _ : state)
        {
            graph.update_routes(elevation);
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
    using rook_cacheall = fs::raster_grid_xt<xt_selector, raster_connect::rook, neighbors_cache<4>>;

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


    BENCHMARK_TEMPLATE(flow_routing__single_flow_router__profile, profile_nocache)
        ->Apply(bms::grid_sizes<benchmark::kMicrosecond>);
    BENCHMARK_TEMPLATE(flow_routing__single_flow_router__profile, profile_cacheall)
        ->Apply(bms::grid_sizes<benchmark::kMicrosecond>);

    BENCH_RASTER(flow_routing__single_flow_router__raster)

}  // namespace fastscapelib::bench
