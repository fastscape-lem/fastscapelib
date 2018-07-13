#include <array>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

#include "fastscapelib/flow_routing.hpp"
#include "fastscapelib/sinks.hpp"
#include "fastscapelib/utils.hpp"

#include "benchmark/benchmark.h"

#include "benchmark_setup.hpp"


namespace bms = benchmark_setup;
namespace fs = fastscapelib;


namespace fastscapelib
{

namespace benchmark_resolve_sinks
{

enum class Method {barnes2014_sloped, kurskal_sloped, boruvka_sloped};


template<bms::SurfaceType surf_type, class T>
class ResolveSinks
{
public:
    ResolveSinks(int n) : s(256) {
        s = bms::FastscapeSetupBase<surf_type, T>(n);
    }

    void resolve_barnes2014_sloped()
    {
        elevation_no_sink = s.elevation;
        fs::fill_sinks_sloped(elevation_no_sink);
        fs::compute_receivers_d8(s.receivers, s.dist2receivers,
                                 elevation_no_sink, s.active_nodes,
                                 s.dx, s.dy);
        fs::compute_donors(s.ndonors, s.donors, s.receivers);
        fs::compute_stack(s.stack, s.ndonors,
                          s.donors, s.receivers);
    }

private:
    bms::FastscapeSetupBase<surf_type, T> s;
    xt::xtensor<T, 2> elevation_no_sink;
};


template<Method method, bms::SurfaceType surf_type, class T>
inline auto bm_resolve_sinks(benchmark::State& state)
{
    auto rs = ResolveSinks<surf_type, T>(state.range(0));

    std::function<void(void)> resolve_func;;

    switch (method) {
    case Method::barnes2014_sloped:
        resolve_func = [&](){ rs.resolve_barnes2014_sloped(); };
        break;
    }

    for (auto _ : state)
    {
        resolve_func();
    }
}


BENCHMARK_TEMPLATE(bm_resolve_sinks,
                   Method::barnes2014_sloped,
                   bms::SurfaceType::cone,
                   double)
->Apply(bms::grid_sizes);

BENCHMARK_TEMPLATE(bm_resolve_sinks,
                   Method::barnes2014_sloped,
                   bms::SurfaceType::cone_noise,
                   double)
->Apply(bms::grid_sizes);

} // namespace benchmark_resolve_sinks

} // namespace fastscapelib
