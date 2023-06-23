#include "benchmark_setup.hpp"

#include "fastscapelib/utils/consts.hpp"
#include "fastscapelib/grid/profile_grid.hpp"
#include "fastscapelib/grid/raster_grid.hpp"

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
        void raster_grid__ctor(benchmark::State& state)
        {
            using grid = G;
            using size_type = typename grid::size_type;

            auto n = static_cast<size_type>(state.range(0));
            std::array<size_type, 2> shape{ { n, n } };

            for (auto _ : state)
            {
                auto rg = grid(shape, { 1., 1. }, fs::node_status::fixed_value);
            }
        }

        template <class G>
        void raster_grid__neighbors_indices(benchmark::State& state)
        {
            using grid_type = G;
            using size_type = typename grid_type::size_type;
            using neighbors_indices_type = typename grid_type::neighbors_indices_type;

            auto n = static_cast<size_type>(state.range(0));
            std::array<size_type, 2> shape{ { n, n } };
            auto grid = grid_type(shape, { 1., 1. }, fs::node_status::fixed_value);

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
        void raster_grid__neighbors_indices_raster(benchmark::State& state)
        {
            using grid_type = G;
            using size_type = typename grid_type::size_type;
            using neighbors_indices_raster_type = typename grid_type::neighbors_indices_raster_type;

            auto n = static_cast<size_type>(state.range(0));
            std::array<size_type, 2> shape{ { n, n } };
            auto grid = grid_type(shape, { 1., 1. }, fs::node_status::fixed_value);

            neighbors_indices_raster_type neighbors_indices;

            // warm-up cache
            for (size_type r = 0; r < shape[0]; ++r)
            {
                for (size_type c = 0; c < shape[1]; ++c)
                {
                    grid.neighbors_indices(r, c, neighbors_indices);
                }
            }

            for (auto _ : state)
            {
                for (size_type r = 0; r < shape[0]; ++r)
                {
                    for (size_type c = 0; c < shape[1]; ++c)
                    {
                        grid.neighbors_indices(r, c, neighbors_indices);
                    }
                }
            }
        }

        template <class G>
        void raster_grid__neighbors(benchmark::State& state)
        {
            using grid_type = G;
            using size_type = typename grid_type::size_type;
            using neighbors_type = typename grid_type::neighbors_type;

            auto n = static_cast<size_type>(state.range(0));
            std::array<size_type, 2> shape{ { n, n } };
            auto grid = grid_type(shape, { 1., 1. }, fs::node_status::fixed_value);

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

        template <class G>
        void raster_grid__neighbors_raster(benchmark::State& state)
        {
            using grid_type = G;
            using size_type = typename grid_type::size_type;
            using neighbors_raster_type = typename grid_type::neighbors_raster_type;

            auto n = static_cast<size_type>(state.range(0));
            std::array<size_type, 2> shape{ { n, n } };
            auto grid = grid_type(shape, { 1., 1. }, fs::node_status::fixed_value);

            neighbors_raster_type neighbors;

            // warm-up cache
            for (size_type r = 0; r < shape[0]; ++r)
            {
                for (size_type c = 0; c < shape[1]; ++c)
                {
                    grid.neighbors(r, c, neighbors);
                }
            }

            for (auto _ : state)
            {
                for (size_type r = 0; r < shape[0]; ++r)
                {
                    for (size_type c = 0; c < shape[1]; ++c)
                    {
                        grid.neighbors(r, c, neighbors);
                    }
                }
            }
        }

        template <class G>
        void raster_grid__neighbor_classic(benchmark::State& state)
        {
            using grid_type = G;
            using size_type = typename grid_type::size_type;
            using neighbors_offsets_type =
                typename raster_neighbors<raster_connect::queen>::neighbors_offsets_type;

            auto n = static_cast<size_type>(state.range(0));
            std::array<size_type, 2> shape{ { n, n } };

            neighbors_offsets_type::shape_type sh0 = { grid_type::n_neighbors_max() };
            neighbors_offsets_type offsets
                = xt::empty<xt::xtensor_fixed<std::ptrdiff_t, xt::xshape<2>>>(sh0);

            auto get_neighbors_indices
                = [&shape, &offsets](auto& r, auto& c) -> neighbors_offsets_type
            {
                for (std::size_t k = 1; k <= grid_type::n_neighbors_max(); ++k)
                {
                    const index_t kr = r + fs::consts::d8_row_offsets[k];
                    const index_t kc = c + fs::consts::d8_col_offsets[k];

                    offsets(k - 1) = xt::xtensor_fixed<std::ptrdiff_t, xt::xshape<2>>({ kr, kc });

                    if (!fs::detail::in_bounds(shape, kr, kc))
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
    BENCHMARK_TEMPLATE(NAME, GRID)->Arg(512)->Unit(benchmark::kMicrosecond);

#define BENCH_RC(NAME)                                                                             \
    BENCH_GRID(NAME, queen_nocache)                                                                \
    BENCH_GRID(NAME, rook_nocache)                                                                 \
    BENCH_GRID(NAME, bishop_nocache)

#define BENCH_ALL(NAME)                                                                            \
    BENCH_GRID(NAME, queen_nocache)                                                                \
    BENCH_GRID(NAME, queen_cacheall)                                                               \
    BENCH_GRID(NAME, rook_nocache)                                                                 \
    BENCH_GRID(NAME, rook_cacheall)                                                                \
    BENCH_GRID(NAME, bishop_nocache)                                                               \
    BENCH_GRID(NAME, bishop_cacheall)


        BENCH_ALL(raster_grid__ctor);
        BENCH_ALL(raster_grid__neighbors_indices);
        BENCH_ALL(raster_grid__neighbors_indices_raster);
        BENCH_ALL(raster_grid__neighbors);
        BENCH_ALL(raster_grid__neighbors_raster);
        BENCH_RC(raster_grid__neighbor_classic);

    }  // namespace bench
}  // namespace fastscapelib
