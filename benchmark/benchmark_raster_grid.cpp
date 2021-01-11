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

        template <class G>
        void raster_grid__ctor(benchmark::State& state)
        {
            using grid = G;
            using size_type = typename grid::size_type;

            auto n = static_cast<size_type>(state.range(0));
            std::array<size_type, 2> shape {n, n};

            for (auto _ : state)
            {
                auto rg = grid(shape, {1., 1.}, fs::node_status::fixed_value_boundary);
            }
        }

        template <class G>
        void raster_grid__neighbors(benchmark::State& state)
        {
            using grid_type = G;
            using size_type = typename grid_type::size_type;
            using neighbors_type = typename grid_type::neighbors_type;

            auto n = static_cast<size_type>(state.range(0));
            std::array<size_type, 2> shape {n, n};
            auto grid = grid_type(shape, {1., 1.}, fs::node_status::fixed_value_boundary);

            // warm-up cache
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

        template <class G>
        void raster_grid__neighbor_view(benchmark::State& state)
        {
            using grid_type = G;
            using size_type = typename grid_type::size_type;

            auto n = static_cast<size_type>(state.range(0));
            std::array<size_type, 2> shape {n, n};

            auto rg = grid_type(shape, {1., 1.}, fs::node_status::fixed_value_boundary);
            xt::xtensor<double, 2> field(shape, 1.);

            for (auto _ : state)
            {
                for (size_type r = 0; r < shape[0]; ++r)
                {
                    for (size_type c = 0; c < shape[1]; ++c)
                    {
                        auto nb_field_view = rg.neighbor_view(field, r * shape[1] + c);
                    }
                }
            }
        }

        template <class G>
        void raster_grid__neighbor_classic(benchmark::State& state)
        {
            using grid_type = G;
            using size_type = typename grid_type::size_type;
            using neighbors_offsets_type = typename raster_neighbors<raster_connect::queen>::neighbors_offsets_type;

            auto n = static_cast<size_type>(state.range(0));
            std::array<size_type, 2> shape {n, n};

            neighbors_offsets_type::shape_type sh0 = {grid_type::max_neighbors()};
            neighbors_offsets_type offsets = xt::empty<xt::xtensor_fixed<std::ptrdiff_t, xt::xshape<2>>>(sh0);

            auto get_neighbors_indices = [&shape, &offsets](auto& r, auto& c) -> neighbors_offsets_type {

                    for(std::size_t k=1; k<=grid_type::max_neighbors(); ++k)
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

        using queen_nocache = fs::raster_grid_xt<xtensor_selector, raster_connect::queen, detail::neighbors_no_cache<8>>;
        using queen_cacheall = fs::raster_grid_xt<xtensor_selector, raster_connect::queen, detail::neighbors_cache<8>>;

        using rook_nocache = fs::raster_grid_xt<xtensor_selector, raster_connect::rook, detail::neighbors_no_cache<4>>;
        using rook_cacheall = fs::raster_grid_xt<xtensor_selector, raster_connect::rook, detail::neighbors_cache<4>>;

        using bishop_nocache = fs::raster_grid_xt<xtensor_selector, raster_connect::bishop, detail::neighbors_no_cache<4>>;
        using bishop_cacheall = fs::raster_grid_xt<xtensor_selector, raster_connect::bishop, detail::neighbors_cache<4>>;


#define BENCH_GRID(NAME, GRID)                              \
    BENCHMARK_TEMPLATE(NAME, GRID)                          \
    ->Apply(bms::small_grid_sizes<benchmark::kMillisecond>);  

#define BENCH_RC(NAME)                \
    BENCH_GRID(NAME, queen_nocache)   \
    BENCH_GRID(NAME, rook_nocache)    \
    BENCH_GRID(NAME, bishop_nocache)

#define BENCH_ALL(NAME)               \
    BENCH_GRID(NAME, queen_nocache)   \
    BENCH_GRID(NAME, queen_cacheall)  \
    BENCH_GRID(NAME, rook_nocache)    \
    BENCH_GRID(NAME, rook_cacheall)   \
    BENCH_GRID(NAME, bishop_nocache)  \
    BENCH_GRID(NAME, bishop_cacheall)


        BENCH_ALL(raster_grid__ctor);
        BENCH_ALL(raster_grid__neighbors);
        BENCH_RC(raster_grid__neighbor_classic);
        BENCH_ALL(raster_grid__neighbor_view);

    } // namespace bench
} // namespace fastscapelib