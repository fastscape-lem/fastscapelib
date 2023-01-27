#include <benchmark/benchmark.h>

#include "fastscapelib/grid/raster_grid.hpp"
#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"

#include "benchmark_setup.hpp"


namespace fs = fastscapelib;
namespace bms = fs::bench_setup;


namespace fastscapelib
{
    namespace bench
    {

        using surf_type = bms::surface_type;


        template <surf_type S, fs::mst_method M>
        void basin_graph__update_routes(benchmark::State& state)
        {
            using flow_graph_type = fs::flow_graph<fs::raster_grid>;
            using basin_graph_type = fs::basin_graph<flow_graph_type::impl_type>;

            auto topo = bms::synthetic_topography_2d<S, double>(state.range(0));
            auto grid = topo.grid();
            auto fgraph = flow_graph_type(grid, { fs::single_flow_router() });

            basin_graph_type bgraph(fgraph.impl(), M);

            auto elevation = topo.elevation();
            fgraph.update_routes(elevation);
            // trigger basin computation
            fgraph.basins();

            for (auto _ : state)
            {
                bgraph.update_routes(elevation);
            }
        }


#define BENCH_SURF_METH(NAME, SURF, METH)                                                          \
    BENCHMARK_TEMPLATE2(NAME, SURF, METH)->Apply(bms::grid_sizes<benchmark::kMillisecond>);

#define BENCH_ALL_METH(NAME, SURF)                                                                 \
    BENCH_SURF_METH(NAME, SURF, fs::mst_method::kruskal)                                           \
    BENCH_SURF_METH(NAME, SURF, fs::mst_method::boruvka)

#define BENCH_ALL_SURF(NAME)                                                                       \
    BENCH_ALL_METH(NAME, surf_type::cone)                                                          \
    BENCH_ALL_METH(NAME, surf_type::cone_inv)                                                      \
    BENCH_ALL_METH(NAME, surf_type::cone_noise)                                                    \
    BENCH_ALL_METH(NAME, surf_type::flat_noise)


        BENCH_ALL_SURF(basin_graph__update_routes);

    }
}
