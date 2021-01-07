#include "benchmark_setup.hpp"

#include "fastscapelib/consts.hpp"
#include "fastscapelib/profile_grid.hpp"
#include "fastscapelib/raster_grid.hpp"

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
        void raster_grid__ctor(benchmark::State& state)
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
        void raster_grid__neighbor_indices(benchmark::State& state)
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

        template <class XT>
        void raster_grid__gcode(benchmark::State& state)
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
                        auto code = rg.gcode(r, c);

                        benchmark::DoNotOptimize(code);
                    }
                }
            }
        }

        template <fs::raster_connect RC>
        void raster_grid__neighbor_offsets(benchmark::State& state)
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
                        auto offset = rg.neighbor_offsets<RC>(0);

                        benchmark::DoNotOptimize(offset);
                    }
                }
            }
        }

        template <fs::raster_connect RC>
        void raster_grid__neighbor_view(benchmark::State& state)
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
        void raster_grid__neighbor_distances(benchmark::State& state)
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
        void raster_grid__neighbor_classic(benchmark::State& state)
        {
            using size_type = typename fs::raster_grid::size_type;

            auto n = static_cast<size_type>(state.range(0));
            std::array<size_type, 2> shape {n, n};

            auto rg = fs::raster_grid(shape, {1., 1.}, fs::node_status::fixed_value_boundary);

            auto get_neighbors_indices = [&shape](auto& r, auto& c) -> fs::raster_grid::offset_list {

                fs::raster_grid::offset_list::shape_type sh0 = {8};
                fs::raster_grid::offset_list offsets = xt::empty<xt::xtensor_fixed<std::ptrdiff_t, xt::xshape<2>>>(sh0);

                for(std::size_t k=1; k<=8; ++k)
                    {
                        const index_t kr = r + fs::consts::d8_row_offsets[k];
                        const index_t kc = c + fs::consts::d8_col_offsets[k];

                        offsets(k-1) = xt::xtensor_fixed<std::ptrdiff_t, xt::xshape<2>>({kr, kc});

                        if(!fs::detail::in_bounds(shape, kr, kc))
                        {
                            continue;
                        }
                    }
                return offsets;
                };

            for (auto _ : state)
            {
                for (size_type r = 0; r < shape[0]; ++r)
                {
                    for (size_type c = 0; c < shape[1]; ++c)
                    {
                        auto res = get_neighbors_indices(r, c);
                    }
                }
            }
        }

        BENCHMARK_TEMPLATE(raster_grid__neighbor_classic, fs::raster_connect::queen)
        ->Apply(bms::grid_sizes<benchmark::kMillisecond>);

        BENCHMARK_TEMPLATE(raster_grid__ctor, fs::xtensor_selector)
        ->Apply(bms::grid_sizes<benchmark::kMicrosecond>);

        BENCHMARK_TEMPLATE(raster_grid__gcode, fs::xtensor_selector)
        ->Apply(bms::grid_sizes<benchmark::kMicrosecond>);

        BENCHMARK_TEMPLATE(raster_grid__neighbor_offsets, fs::raster_connect::queen)
        ->Apply(bms::grid_sizes<benchmark::kMillisecond>);

        BENCHMARK_TEMPLATE(raster_grid__neighbor_indices, fs::raster_connect::queen)
        ->Apply(bms::grid_sizes<benchmark::kMillisecond>);

        BENCHMARK_TEMPLATE(raster_grid__neighbor_indices, fs::raster_connect::rook)
        ->Apply(bms::grid_sizes<benchmark::kMillisecond>);

        BENCHMARK_TEMPLATE(raster_grid__neighbor_view, fs::raster_connect::queen)
        ->Apply(bms::grid_sizes<benchmark::kMillisecond>);

        BENCHMARK_TEMPLATE(raster_grid__neighbor_view, fs::raster_connect::rook)
        ->Apply(bms::grid_sizes<benchmark::kMillisecond>);

        BENCHMARK_TEMPLATE(raster_grid__neighbor_distances, fs::raster_connect::queen)
        ->Apply(bms::grid_sizes<benchmark::kMillisecond>);

        BENCHMARK_TEMPLATE(raster_grid__neighbor_distances, fs::raster_connect::rook)
        ->Apply(bms::grid_sizes<benchmark::kMillisecond>);

    } // namespace bench
} // namespace fastscapelib