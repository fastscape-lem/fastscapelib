#include "benchmark_setup.hpp"

#include "fastscapelib/consts.hpp"
#include "fastscapelib/profile_grid.hpp"

#include "xtensor/xfixed.hpp"
#include "xtensor/xtensor.hpp"

#include <benchmark/benchmark.h>


namespace fs = fastscapelib;
namespace bms = benchmark_setup;


namespace fastscapelib
{
    namespace bench
    {

        template <class XT>
        void profile_grid__ctor(benchmark::State& state)
        {
            using grid_type = fs::profile_grid_xt<XT>;
            using size_type = typename grid_type::size_type;

            auto shape = static_cast<size_type>(state.range(0));

            for (auto _ : state)
            {
                auto grid = grid_type(shape, 1.3, fs::node_status::fixed_value_boundary);
            }
        }

        template <class XT>
        void profile_grid__neighbors(benchmark::State& state)
        {
            using grid_type = fs::profile_grid_xt<XT>;
            using size_type = typename grid_type::size_type;
            using neighbors_type = typename grid_type::neighbors_type;

            auto size = static_cast<size_type>(state.range(0));
            auto grid = grid_type(size, 1.3, fs::node_status::fixed_value_boundary);

            //warm-up cache
            for (size_type idx = 0; idx < grid.size(); ++idx)
            {
                auto neighbors = grid.neighbors(idx);
            } 

            neighbors_type neighbors;
            
            for (auto _ : state)
            {
                for (size_type idx = 0; idx < grid.size(); ++idx)
                {
                    grid.neighbors(idx, neighbors);
                }
            }
        }

        template <class XT>
        void profile_grid__spacing(benchmark::State& state)
        {
            using grid_type = fs::profile_grid_xt<XT>;
            using size_type = typename grid_type::size_type;

            auto size = static_cast<size_type>(state.range(0));
            auto grid = grid_type(size, 1.3, fs::node_status::fixed_value_boundary);

            for (auto _ : state)
            {
                auto spacing = grid.spacing();
                benchmark::DoNotOptimize(spacing);
            }
        }

        template <class XT>
        void profile_grid__size(benchmark::State& state)
        {
            using grid_type = fs::profile_grid_xt<XT>;
            using size_type = typename grid_type::size_type;

            auto size = static_cast<size_type>(state.range(0));
            auto grid = grid_type(size, 1.3, fs::node_status::fixed_value_boundary);

            for (auto _ : state)
            {
                auto s = grid.size();
                benchmark::DoNotOptimize(s);
            }
        }

        BENCHMARK_TEMPLATE(profile_grid__ctor, fs::xtensor_selector)
        ->Apply(bms::grid_sizes<benchmark::kMicrosecond>);

        BENCHMARK_TEMPLATE(profile_grid__neighbors, fs::xtensor_selector)
        ->Apply(bms::grid_sizes<benchmark::kMicrosecond>);

        BENCHMARK_TEMPLATE(profile_grid__spacing, fs::xtensor_selector)
        ->Apply(bms::grid_sizes<benchmark::kNanosecond>);

        BENCHMARK_TEMPLATE(profile_grid__size, fs::xtensor_selector)
        ->Apply(bms::grid_sizes<benchmark::kNanosecond>);
 
    } // namespace bench
} // namespace fastscapelib