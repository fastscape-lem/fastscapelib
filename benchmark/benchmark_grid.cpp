#include <benchmark/benchmark.h>

#include "xtensor/xfixed.hpp"
#include "xtensor/xtensor.hpp"

#include "fastscapelib/consts.hpp"
#include "fastscapelib/grid.hpp"

#include "benchmark_setup.hpp"


namespace fs = fastscapelib;
namespace bms = benchmark_setup;


namespace fastscapelib
{

namespace benchmark_grid
{

template <class XT>
void bm_raster_grid_create(benchmark::State& state)
{
    using grid = fs::raster_grid_xt<XT>;
    using size_type = typename grid::size_type;

    auto n = static_cast<size_type>(state.range(0));
    std::array<size_type, 2> shape {n, n};

    for (auto _ : state)
    {
        auto rg = grid(shape, {1., 1.}, fs::node_status::fixed_value_boundary);
    }
}

template <fs::raster_connect RC>
void bm_raster_grid_neighbor_indices(benchmark::State& state)
{
    using size_type = typename fs::raster_grid::size_type;

    auto n = static_cast<size_type>(state.range(0));
    std::array<size_type, 2> shape {n, n};

    auto rg = fs::raster_grid(shape, {1., 1.}, fs::node_status::fixed_value_boundary);

    for (auto _ : state)
    {
        for (size_type r = 0; r < shape[0]; ++r)
        {
            for (size_type c = 0; c < shape[1]; ++c)
            {
                auto nb_indices = rg.neighbor_indices<RC>(r, c);
            }
        }
    }
}

template <fs::raster_connect RC>
void bm_raster_grid_neighbor_view(benchmark::State& state)
{
    using size_type = typename fs::raster_grid::size_type;

    auto n = static_cast<size_type>(state.range(0));
    std::array<size_type, 2> shape {n, n};

    auto rg = fs::raster_grid(shape, {1., 1.}, fs::node_status::fixed_value_boundary);
    xt::xtensor<double, 2> field(shape, 1.);

    for (auto _ : state)
    {
        for (size_type r = 0; r < shape[0]; ++r)
        {
            for (size_type c = 0; c < shape[1]; ++c)
            {
                auto nb_field_view = rg.neighbor_view<RC>(field, r, c);
            }
        }
    }
}

template <fs::raster_connect RC>
void bm_raster_grid_neighbor_distances(benchmark::State& state)
{
    using size_type = typename fs::raster_grid::size_type;

    auto n = static_cast<size_type>(state.range(0));
    std::array<size_type, 2> shape {n, n};

    auto rg = fs::raster_grid(shape, {1., 1.}, fs::node_status::fixed_value_boundary);

    for (auto _ : state)
    {
        for (size_type r = 0; r < shape[0]; ++r)
        {
            for (size_type c = 0; c < shape[1]; ++c)
            {
                auto nb_distances = rg.neighbor_distances<RC>(r, c);
            }
        }
    }
}

template <fs::raster_connect RC>
void bm_raster_grid_neighbor_classic(benchmark::State& state)
{
    using size_type = typename fs::raster_grid::size_type;

    auto n = static_cast<size_type>(state.range(0));
    std::array<size_type, 2> shape {n, n};

    auto rg = fs::raster_grid(shape, {1., 1.}, fs::node_status::fixed_value_boundary);

    for (auto _ : state)
    {
        for (size_type r = 0; r < shape[0]; ++r)
        {
            for (size_type c = 0; c < shape[1]; ++c)
            {
                for(std::size_t k=1; k<=8; ++k)
                {
                    const index_t kr = r + fs::consts::d8_row_offsets[k];
                    const index_t kc = c + fs::consts::d8_col_offsets[k];

                    if(!fs::detail::in_bounds(shape, kr, kc))
                    {
                        continue;
                    }

                    benchmark::DoNotOptimize(kr);
                    benchmark::DoNotOptimize(kc);
                }
            }
        }
    }
}

BENCHMARK_TEMPLATE(bm_raster_grid_neighbor_classic, fs::raster_connect::queen)
->Apply(bms::grid_sizes<benchmark::kMillisecond>);

BENCHMARK_TEMPLATE(bm_raster_grid_create, fs::xtensor_selector)
->Apply(bms::grid_sizes<benchmark::kMicrosecond>);

BENCHMARK_TEMPLATE(bm_raster_grid_neighbor_indices, fs::raster_connect::queen)
->Apply(bms::grid_sizes<benchmark::kMillisecond>);

BENCHMARK_TEMPLATE(bm_raster_grid_neighbor_indices, fs::raster_connect::rook)
->Apply(bms::grid_sizes<benchmark::kMillisecond>);

BENCHMARK_TEMPLATE(bm_raster_grid_neighbor_view, fs::raster_connect::queen)
->Apply(bms::grid_sizes<benchmark::kMillisecond>);

BENCHMARK_TEMPLATE(bm_raster_grid_neighbor_view, fs::raster_connect::rook)
->Apply(bms::grid_sizes<benchmark::kMillisecond>);

BENCHMARK_TEMPLATE(bm_raster_grid_neighbor_distances, fs::raster_connect::queen)
->Apply(bms::grid_sizes<benchmark::kMillisecond>);

BENCHMARK_TEMPLATE(bm_raster_grid_neighbor_distances, fs::raster_connect::rook)
->Apply(bms::grid_sizes<benchmark::kMillisecond>);

} // namespace benchmark_grid

} // namespace fastscapelib
