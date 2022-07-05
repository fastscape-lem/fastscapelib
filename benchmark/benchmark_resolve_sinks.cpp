#include <array>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xnpy.hpp"

#include "fastscapelib/flow_routing.hpp"
#include "fastscapelib/sinks.hpp"
#include "fastscapelib/utils.hpp"

#ifdef ENABLE_RICHDEM
#include "fastscapelib/richdem.hpp"
#endif

#include "benchmark/benchmark.h"

#include "benchmark_setup.hpp"


namespace bms = benchmark_setup;
namespace fs = fastscapelib;

using SurfaceType = bms::SurfaceType;


namespace fastscapelib
{

    namespace benchmark_resolve_sinks
    {

        enum class Method
        {
            barnes2014_sloped,
            kruskal_sloped,
            boruvka_sloped,
            flats_then_sloped
        };


        template <SurfaceType surf_type, class T>
        class ResolveSinks
        {
        public:
            ResolveSinks(int n)
                : s(256)
            {
                s = bms::FastscapeSetupBase<surf_type, T>(n);

                if (surf_type == SurfaceType::custom)
                {
                    s.elevation = xt::load_npy<T>("bedrock_4096.npy");
                }
            }

            void resolve_barnes2014_sloped()
            {
                elevation_no_sink = s.elevation;
                fs::fill_sinks_sloped(elevation_no_sink);

                fs::compute_receivers_d8(
                    s.receivers, s.dist2receivers, elevation_no_sink, s.active_nodes, s.dx, s.dy);
                fs::compute_donors(s.ndonors, s.donors, s.receivers);
                fs::compute_stack(s.stack, s.ndonors, s.donors, s.receivers);
            }

            template <fs::mst_method basin_algo, fs::sink_route_method connect_type>
            void resolve_basin_graph()
            {
                fs::compute_receivers_d8(
                    s.receivers, s.dist2receivers, s.elevation, s.active_nodes, s.dx, s.dy);
                fs::compute_donors(s.ndonors, s.donors, s.receivers);
                fs::compute_stack(s.stack, s.ndonors, s.donors, s.receivers);

                fs::correct_flowrouting<basin_algo, connect_type>(basin_graph,
                                                                  s.basins,
                                                                  s.receivers,
                                                                  s.dist2receivers,
                                                                  s.ndonors,
                                                                  s.donors,
                                                                  s.stack,
                                                                  s.active_nodes,
                                                                  s.elevation,
                                                                  s.dx,
                                                                  s.dy);
            }

#ifdef ENABLE_RICHDEM
            void resolve_flats_then_sloped()
            {
                elevation_no_sink = s.elevation;
                fs::fill_sinks_wei2018(elevation_no_sink);
                fs::resolve_flats_sloped(elevation_no_sink);

                fs::compute_receivers_d8(
                    s.receivers, s.dist2receivers, elevation_no_sink, s.active_nodes, s.dx, s.dy);
                fs::compute_donors(s.ndonors, s.donors, s.receivers);
                fs::compute_stack(s.stack, s.ndonors, s.donors, s.receivers);
            }
#endif

        private:
            bms::FastscapeSetupBase<surf_type, T> s;
            xt::xtensor<T, 2> elevation_no_sink;
            fastscapelib::BasinGraph<index_t, index_t, T> basin_graph;
        };


        template <Method method, SurfaceType surf_type, class T>
        inline auto bm_resolve_sinks(benchmark::State& state)
        {
            auto rs = ResolveSinks<surf_type, T>(state.range(0));

            std::function<void(void)> resolve_func;
            ;

            switch (method)
            {
                case Method::barnes2014_sloped:
                    resolve_func = [&]() { rs.resolve_barnes2014_sloped(); };
                    break;

                case Method::kruskal_sloped:
                    resolve_func = [&]() {
                        rs.template resolve_basin_graph<fs::mst_method::kruskal,
                                                        fs::sink_route_method::fill_sloped>();
                    };
                    break;

                case Method::boruvka_sloped:
                    resolve_func = [&]() {
                        rs.template resolve_basin_graph<fs::mst_method::boruvka,
                                                        fs::sink_route_method::fill_sloped>();
                    };
                    break;

                case Method::flats_then_sloped:
                    resolve_func = [&]() { rs.resolve_flats_then_sloped(); };
                    break;
            }

            for (auto _ : state)
            {
                resolve_func();
            }
        }


        BENCHMARK_TEMPLATE(bm_resolve_sinks, Method::barnes2014_sloped, SurfaceType::custom, double)
            ->Args({ 4096 })
            ->Unit(benchmark::kMillisecond);

        BENCHMARK_TEMPLATE(bm_resolve_sinks, Method::kruskal_sloped, SurfaceType::custom, double)
            ->Args({ 4096 })
            ->Unit(benchmark::kMillisecond);

        BENCHMARK_TEMPLATE(bm_resolve_sinks, Method::boruvka_sloped, SurfaceType::custom, double)
            ->Args({ 4096 })
            ->Unit(benchmark::kMillisecond);

        BENCHMARK_TEMPLATE(bm_resolve_sinks, Method::flats_then_sloped, SurfaceType::custom, double)
            ->Args({ 4096 })
            ->Unit(benchmark::kMillisecond);

        BENCHMARK_TEMPLATE(bm_resolve_sinks, Method::barnes2014_sloped, SurfaceType::cone, double)
            ->Apply(bms::grid_sizes);

        BENCHMARK_TEMPLATE(bm_resolve_sinks,
                           Method::barnes2014_sloped,
                           SurfaceType::cone_noise,
                           double)
            ->Apply(bms::grid_sizes);

        BENCHMARK_TEMPLATE(bm_resolve_sinks,
                           Method::barnes2014_sloped,
                           SurfaceType::flat_noise,
                           double)
            ->Apply(bms::grid_sizes);

        BENCHMARK_TEMPLATE(bm_resolve_sinks, Method::kruskal_sloped, SurfaceType::cone, double)
            ->Apply(bms::grid_sizes);

        BENCHMARK_TEMPLATE(bm_resolve_sinks,
                           Method::kruskal_sloped,
                           SurfaceType::cone_noise,
                           double)
            ->Apply(bms::grid_sizes);


        BENCHMARK_TEMPLATE(bm_resolve_sinks,
                           Method::kruskal_sloped,
                           SurfaceType::flat_noise,
                           double)
            ->Apply(bms::grid_sizes);

        BENCHMARK_TEMPLATE(bm_resolve_sinks, Method::boruvka_sloped, SurfaceType::cone, double)
            ->Apply(bms::grid_sizes);

        BENCHMARK_TEMPLATE(bm_resolve_sinks,
                           Method::boruvka_sloped,
                           SurfaceType::cone_noise,
                           double)
            ->Apply(bms::grid_sizes);


        BENCHMARK_TEMPLATE(bm_resolve_sinks,
                           Method::boruvka_sloped,
                           SurfaceType::flat_noise,
                           double)
            ->Apply(bms::grid_sizes);

        BENCHMARK_TEMPLATE(bm_resolve_sinks, Method::flats_then_sloped, SurfaceType::cone, double)
            ->Apply(bms::grid_sizes);

        BENCHMARK_TEMPLATE(bm_resolve_sinks,
                           Method::flats_then_sloped,
                           SurfaceType::cone_noise,
                           double)
            ->Apply(bms::grid_sizes);

        BENCHMARK_TEMPLATE(bm_resolve_sinks,
                           Method::flats_then_sloped,
                           SurfaceType::flat_noise,
                           double)
            ->Apply(bms::grid_sizes);


    }  // namespace benchmark_resolve_sinks

}  // namespace fastscapelib
