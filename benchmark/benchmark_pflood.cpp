#include <benchmark/benchmark.h>

#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"

#include "benchmark_setup.hpp"


namespace fs = fastscapelib;
namespace bms = fs::bench_setup;


namespace fastscapelib
{
    class dummy_flow_router : public flow_operator
    {
    public:
        inline std::string name() const noexcept override
        {
            return "dummy_flow_router";
        }

        static constexpr bool graph_updated = true;
        static constexpr flow_direction out_flowdir = flow_direction::single;
    };


    namespace detail
    {
        template <class FG>
        class flow_operator_impl<FG, dummy_flow_router, flow_graph_fixed_array_tag>
            : public flow_operator_impl_base<FG, dummy_flow_router>
        {
        public:
            using graph_impl_type = FG;
            using base_type = flow_operator_impl_base<FG, dummy_flow_router>;
            using data_array_type = typename graph_impl_type::data_array_type;

            flow_operator_impl(std::shared_ptr<dummy_flow_router> ptr)
                : base_type(std::move(ptr)){};

            void apply(graph_impl_type& /*graph_impl*/, data_array_type& /*elevation*/)
            {
            }
        };
    }
}

namespace fastscapelib
{
    namespace bench
    {

        using surf_type = bms::surface_type;

        template <surf_type S>
        void fill_sinks_sloped(benchmark::State& state)
        {
            using topography_type = typename bms::synthetic_topography_2d<S, double>;
            using grid_type = typename topography_type::grid_type;
            using size_type = typename grid_type::size_type;

            topography_type topo(state.range(0));
            auto grid = topo.grid();

            xt::xtensor<double, 2> elevation = xt::random::rand<double>(grid.shape());
            flow_graph<grid_type> graph(grid, { pflood_sink_resolver(), dummy_flow_router() });

            // warm-up cache
            for (size_type idx = 0; idx < grid.size(); ++idx)
            {
                grid.neighbors(idx);
            }

            for (auto _ : state)
            {
                graph.update_routes(elevation);
            }
        }


#define BENCH_SURF(NAME, SURF)                                                                     \
    BENCHMARK_TEMPLATE(NAME, SURF)->Apply(bms::small_grid_sizes<benchmark::kMillisecond>);

#define BENCH_ALL_SURF(NAME)                                                                       \
    BENCH_SURF(NAME, surf_type::cone)                                                              \
    BENCH_SURF(NAME, surf_type::cone_inv)                                                          \
    BENCH_SURF(NAME, surf_type::cone_noise)                                                        \
    BENCH_SURF(NAME, surf_type::flat_noise)


        BENCH_ALL_SURF(fill_sinks_sloped);

    }
}
