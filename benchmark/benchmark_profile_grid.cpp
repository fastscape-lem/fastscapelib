#include "benchmark_setup.hpp"

#include "fastscapelib/utils/consts.hpp"
#include "fastscapelib/grid/profile_grid.hpp"

#include "xtensor/xfixed.hpp"
#include "xtensor/xtensor.hpp"

#include <benchmark/benchmark.h>


namespace fs = fastscapelib;
namespace bms = fs::bench_setup;


namespace fastscapelib
{
    namespace bench
    {

        template <class G>
        void profile_grid__ctor(benchmark::State& state)
        {
            using grid_type = G;
            using size_type = typename grid_type::size_type;

            auto shape = static_cast<size_type>(state.range(0));

            for (auto _ : state)
            {
                auto grid = grid_type(shape, 1.3, fs::node_status::fixed_value);
            }
        }

        template <class G>
        void profile_grid__neighbors_indices(benchmark::State& state)
        {
            using grid_type = G;
            using size_type = typename grid_type::size_type;
            using neighbors_indices_type = typename grid_type::neighbors_indices_type;

            auto size = static_cast<size_type>(state.range(0));
            auto grid = grid_type(size, 1.3, fs::node_status::fixed_value);

            neighbors_indices_type neighbors_indices;

            // warm-up cache
            for (size_type idx = 0; idx < grid.size(); ++idx)
            {
                grid.neighbors_indices(idx, neighbors_indices);
            }

            for (auto _ : state)
            {
                for (size_type idx = 0; idx < grid.size(); ++idx)
                {
                    grid.neighbors_indices(idx, neighbors_indices);
                }
            }
        }

        template <class G>
        void profile_grid__neighbors(benchmark::State& state)
        {
            using grid_type = G;
            using size_type = typename grid_type::size_type;
            using neighbors_type = typename grid_type::neighbors_type;

            auto size = static_cast<size_type>(state.range(0));
            auto grid = grid_type(size, 1.3, fs::node_status::fixed_value);

            neighbors_type neighbors;

            // warm-up cache
            for (size_type idx = 0; idx < grid.size(); ++idx)
            {
                grid.neighbors(idx, neighbors);
            }

            for (auto _ : state)
            {
                for (size_type idx = 0; idx < grid.size(); ++idx)
                {
                    grid.neighbors(idx, neighbors);
                }
            }
        }

        using profile_nocache = fs::profile_grid_xt<xt_selector, neighbors_no_cache<2>>;
        using profile_cacheall = fs::profile_grid;

#define BENCH_GRID(NAME, GRID)                                                                     \
    BENCHMARK_TEMPLATE(NAME, GRID)->Arg(512)->Unit(benchmark::kNanosecond);

#define BENCH_RASTER(NAME)                                                                         \
    BENCH_GRID(NAME, profile_nocache)                                                              \
    BENCH_GRID(NAME, profile_cacheall)

        BENCH_RASTER(profile_grid__ctor);
        BENCH_RASTER(profile_grid__neighbors_indices);
        BENCH_RASTER(profile_grid__neighbors);

    }  // namespace bench
}  // namespace fastscapelib
