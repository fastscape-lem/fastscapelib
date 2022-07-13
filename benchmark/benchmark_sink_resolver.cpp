/*
 * Benchmak route flow and resolve solve closed depressions at the same time.
 */
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

        /*
         * Lots of small closed depressions
         */
        template <class SR>
        void sink_resolver__flat_noise(benchmark::State& state, SR&& resolver)
        {
            auto topo = bms::synthetic_topography_2d<surf_type::flat_noise, double>(state.range(0));
            auto grid = topo.grid();
            auto fgraph = fs::make_flow_graph(grid, fs::single_flow_router(), resolver);

            auto elevation = topo.elevation();

            for (auto _ : state)
            {
                fgraph.update_routes(elevation);
            }
        }


        /*
         * No closed depression
         */
        template <class SR>
        void sink_resolver__cone(benchmark::State& state, SR&& resolver)
        {
            auto topo = bms::synthetic_topography_2d<surf_type::cone, double>(state.range(0));
            auto grid = topo.grid();
            auto fgraph = fs::make_flow_graph(grid, fs::single_flow_router(), resolver);

            auto elevation = topo.elevation();

            for (auto _ : state)
            {
                fgraph.update_routes(elevation);
            }
        }


        /*
         * One large closed depression
         */
        template <class SR>
        void sink_resolver__cone_inv(benchmark::State& state, SR&& resolver)
        {
            auto topo = bms::synthetic_topography_2d<surf_type::cone_inv, double>(state.range(0));
            auto grid = topo.grid();
            auto fgraph = fs::make_flow_graph(grid, fs::single_flow_router(), resolver);

            auto elevation = topo.elevation();

            for (auto _ : state)
            {
                fgraph.update_routes(elevation);
            }
        }


#define BENCH_RESOLVER(FUNC, NAME, RESOLVER)                                                       \
    BENCHMARK_CAPTURE(FUNC, NAME, RESOLVER)->Apply(bms::small_grid_sizes<benchmark::kMillisecond>);


#define BENCH_SURF(FUNC, NAME)                                                                     \
    BENCH_RESOLVER(FUNC, NAME - pflood, fs::pflood_sink_resolver())                                \
    BENCH_RESOLVER(                                                                                \
        FUNC,                                                                                      \
        NAME - mst - kruskal - basic,                                                              \
        (fs::mst_sink_resolver{ fs::mst_method::kruskal, fs::mst_route_method::basic }))           \
    BENCH_RESOLVER(                                                                                \
        FUNC,                                                                                      \
        NAME - mst - kruskal - carve,                                                              \
        (fs::mst_sink_resolver{ fs::mst_method::kruskal, fs::mst_route_method::carve }))           \
    BENCH_RESOLVER(                                                                                \
        FUNC,                                                                                      \
        NAME - mst - boruvka - basic,                                                              \
        (fs::mst_sink_resolver{ fs::mst_method::boruvka, fs::mst_route_method::basic }))           \
    BENCH_RESOLVER(                                                                                \
        FUNC,                                                                                      \
        NAME - mst - boruvka - carve,                                                              \
        (fs::mst_sink_resolver{ fs::mst_method::boruvka, fs::mst_route_method::carve }))


        BENCH_SURF(sink_resolver__flat_noise, flat_noise);
        BENCH_SURF(sink_resolver__cone, cone);
        BENCH_SURF(sink_resolver__cone_inv, cone_inv);
    }
}
