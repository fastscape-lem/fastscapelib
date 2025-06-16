#include "benchmark_setup.hpp"

#include "fastscapelib/grid/raster_grid.hpp"
#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"
#include "fastscapelib/utils/eigen_containers.hpp"

#include "xtensor/containers/xtensor.hpp"
#include "xtensor/generators/xrandom.hpp"

#include <benchmark/benchmark.h>


namespace fs = fastscapelib;
namespace bms = fs::bench_setup;


template <class T>
struct random_elevation;


template <class T, std::size_t N>
struct random_elevation<xt::xtensor<T, N>>
{
    static auto generate(std::size_t n)
    {
        return xt::random::rand<double>({ n, n });
    }
};

template <class T, int R, int C>
struct random_elevation<Eigen::Array<T, R, C>>
{
    static auto generate(std::size_t n)
    {
        return Eigen::Array<T, R, C>::Random(n, n);
    }
};

namespace fastscapelib
{
    namespace bench
    {

        template <class G>
        void single_flow_router__raster(benchmark::State& state)
        {
            using grid_type = G;
            using size_type = typename grid_type::size_type;

            auto n = static_cast<size_type>(state.range(0));
            std::array<size_type, 2> shape{ { n, n } };
            auto grid = grid_type(shape, { 1., 1. }, fs::node_status::fixed_value);
            auto graph = fs::flow_graph<grid_type>(grid, { fs::single_flow_router() });

            auto elevation = random_elevation<typename grid_type::container_type>::generate(n);

            // warm-up grid cache
            for (size_type idx = 0; idx < grid.size(); ++idx)
            {
                grid.neighbors(idx);
            }

            for (auto _ : state)
            {
                graph.update_routes(elevation);
            }
        }

        using xt_queen_cacheall
            = fs::raster_grid<xt_selector, raster_connect::queen, neighbors_cache<8>>;
        using eigen_queen_cacheall
            = fs::raster_grid<eigen_selector, raster_connect::queen, neighbors_cache<8>>;

        BENCHMARK_TEMPLATE(single_flow_router__raster, xt_queen_cacheall)
            ->Apply(bms::grid_sizes<benchmark::kMillisecond>);
        // BENCHMARK_TEMPLATE(single_flow_router__raster, eigen_queen_cacheall)
        //     ->Apply(bms::grid_sizes<benchmark::kMillisecond>);

    }  // namespace bench
}  // namespace fastscapelib
