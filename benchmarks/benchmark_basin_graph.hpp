#pragma once

#include "fastscapelib/basin_graph.hpp"
#include "fastscapelib/flow_routing.hpp"

#include "bm_types.hpp"

extern int boruvka_perf;

template <fastscapelib::BasinAlgo algo, fastscapelib::ConnectType connect>
void
benchmark_fastscape_basin(xt::xtensor<bm_index, 1>& stack,
                          xt::xtensor<bm_index, 1>& receivers,
                          xt::xtensor<bm_scalar, 1>& dist2receivers,
                          const xt::xtensor<bm_scalar, 2>& elevation,
                          const xt::xtensor<bool, 2>& active_nodes,
                          bm_scalar dx,
                          bm_scalar dy)
{
    const auto elev_shape = elevation.shape();
    const size_t nrows = (size_t) elev_shape[0];
    const size_t ncols = (size_t) elev_shape[1];
    std::array<std::size_t, 1> shape_1D{ nrows * ncols };

    fastscapelib::BasinGraph<bm_index, bm_index, bm_scalar> basin_graph;

    xt::xtensor<bm_index, 1> basins(shape_1D);
    xt::xtensor<bm_index, 1> ndonors(shape_1D);
    xt::xtensor<bm_index, 2> donors({ nrows * ncols, 8 });


    fastscapelib::compute_receivers_d8(receivers, dist2receivers, elevation, active_nodes, dx, dy);

    fastscapelib::compute_donors(ndonors, donors, receivers);

    fastscapelib::compute_stack(stack, ndonors, donors, receivers);

    boruvka_perf = -1;
    fastscapelib::correct_flowrouting<algo, connect>(basin_graph,
                                                     basins,
                                                     receivers,
                                                     dist2receivers,
                                                     ndonors,
                                                     donors,
                                                     stack,
                                                     active_nodes,
                                                     elevation,
                                                     dx,
                                                     dy);
}
