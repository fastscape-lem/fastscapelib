/**
 * @file
 * @brief Provides implementation of various efficient
 *        algorithms related to flow routing.
 */
#pragma once

#include <cmath>
#include <algorithm>
#include <limits>
#include <array>

#include "xtensor/xtensor.hpp"
#include "xtensor/xadapt.hpp"

#include "fastscapelib/utils.hpp"
#include "fastscapelib/consts.hpp"

#include "fastscapelib/basin_graph.hpp"

#include "fastscapelib/Profile.h"
#include <assert.h>
#include <xtensor/xio.hpp>

namespace fs = fastscapelib;


namespace fastscapelib
{

namespace detail
{

inline auto get_d8_distances_inv(double dx, double dy) -> std::array<double, 9>
{
    std::array<double, 9> d8_dists;

    for(size_t k=0; k<9; ++k)
    {
        double d8_dx = dx * fs::consts::d8_col_offsets[k];
        double d8_dy = dy * fs::consts::d8_row_offsets[k];
        d8_dists[k] = 1.0/std::sqrt(d8_dy*d8_dy + d8_dx*d8_dx);
    }

    return d8_dists;
}


template<class A1, class A2, class A3>
void add2stack(index_t& nstack,
               A1& stack,
               const A2& ndonors,
               const A3& donors,
               index_t inode)
{
    for(unsigned short k=0; k<ndonors(inode); ++k)
    {
        index_t idonor = donors(inode, k);
        stack(nstack) = idonor;
        ++nstack;
        add2stack(nstack, stack, ndonors, donors, idonor);
    }
}

}  // namespace detail

template<class A1, class A2, class A3, class A4>
void compute_receivers_d8(A1& receivers,
                          A2& dist2receivers,
                          const A3& elevation,
                          const A4& active_nodes,
                          double dx,
                          double dy)
{
    using elev_t = typename A3::value_type;
    const auto d8_dists = detail::get_d8_distances_inv(dx, dy);

    const auto elev_shape = elevation.shape();
    const index_t nrows = (index_t) elev_shape[0];
    const index_t ncols = (index_t) elev_shape[1];

    for(index_t r=0; r<nrows; ++r)
    {
        for(index_t c=0; c<ncols; ++c)
        {
            index_t inode = r * ncols + c;

            receivers(inode) = inode;
            dist2receivers(inode) = 0.;

            if(!active_nodes(inode))
                continue;

            double slope_max = 0.0;
            size_t k_found = 0;

            for(size_t k=1; k<=8; ++k)
            {
                index_t kr = r + fs::consts::d8_row_offsets[k];
                index_t kc = c + fs::consts::d8_col_offsets[k];

                if(!fs::detail::in_bounds(elev_shape, kr, kc))
                    continue;

                index_t ineighbor = kr * ncols + kc;
                double slope = (elevation(inode) - elevation(ineighbor)) * d8_dists[k];

                if(slope > slope_max)
                {
                    slope_max = slope;
                    k_found = k;
                }
            }

            index_t kr = r + fs::consts::d8_row_offsets[k_found];
            index_t kc = c + fs::consts::d8_col_offsets[k_found];

            receivers(inode) = kr * ncols + kc;
            dist2receivers(inode) = d8_dists[k_found];
        }
    }
}


template<class A1, class A2, class A3>
void compute_donors(A1& ndonors, A2& donors, const A3& receivers)
{
    index_t nnodes = (index_t) receivers.size();

    std::fill(ndonors.begin(), ndonors.end(), 0);

    for(index_t inode=0; inode<nnodes; ++inode)
    {
        if(receivers(inode) != inode)
        {
            index_t irec = receivers(inode);
            donors(irec, ndonors(irec)) = inode;
            ++ndonors(irec);
        }
    }
}


template<class A1, class A2, class A3, class A4>
void compute_stack(A1& stack,
                   const A2& ndonors,
                   const A3& donors,
                   const A4& receivers)
{
    index_t nnodes = (index_t) receivers.size();
    index_t nstack = 0;

    for(index_t inode=0; inode<nnodes; ++inode)
    {
        if(receivers(inode) == inode)
        {
            stack(nstack) = inode;
            ++nstack;
            detail::add2stack(nstack, stack, ndonors, donors, inode);
        }
    }
}


template <class Basins_XT, class Stack_XT, class Rcv_XT>
auto compute_basins(Basins_XT& basins,
                    const Stack_XT& stack,
                    const Rcv_XT& receivers)
-> typename Basins_XT::value_type         // returns last basin +1
{

    using Basin_T = typename Basins_XT::value_type;
    using Node_T  = typename Stack_XT::value_type;
    BasinGraph<Basin_T, Node_T> basin_graph;

    basin_graph.compute_basins(basins, stack, receivers);
    return basin_graph.basin_count();

}


template<class A1, class A2, class A3>
index_t compute_pits(A1& pits,
                     const A2& outlets,
                     const A3& active_nodes,
                     index_t nbasins)
{
    //TODO: works whether active_nodes is 1-d or 2-d but nif(!active_nodes(r, c)) { continue; }ot safe!
    //      see xtensor issue #588
    index_t ipit = 0;

    for(index_t ibasin=0; ibasin<nbasins; ++ibasin)
    {
        index_t inode = outlets(ibasin);

        if(active_nodes(inode))
        {
            pits(ipit) = inode;
            ++ipit;
        }
    }

    index_t npits = ipit;

    return npits;
}


template<class A1, class A2, class A3>
void compute_drainage_area(A1& area, const A2& stack, const A3& receivers)
{
    for(auto inode=stack.crbegin(); inode!=stack.crend(); ++inode)
    {
        if(receivers(*inode) != *inode)
        {
            area(receivers(*inode)) += area(*inode);
        }
    }
}


template<class A1, class A2, class A3>
void compute_drainage_area(A1& area,
                           const A2& stack,
                           const A3& receivers,
                           double dx,
                           double dy)
{
    std::fill(area.begin(), area.end(), dx * dy);

    //TODO: replace with safe/clean way to get flatten view in xtensor.
    //      (see xtensor issues #322 #324 #588).
    auto area_flat = xt::adapt(area.data(),
                               std::array<size_t, 1>{ stack.size() });

    compute_drainage_area(area_flat, stack, receivers);
}


template <BasinAlgo algo, ConnectType connect, class BasinGraph_T, class Basins_XT, class Rcv_XT,
          class DistRcv_XT, class NDonnors_XT, class Donnors_XT,
          class Stack_XT, class Active_XT,
          class Elevation_XT>
void correct_flowrouting(BasinGraph_T& basin_graph, Basins_XT& basins,
                         Rcv_XT& receivers, DistRcv_XT& dist2receivers,
                         NDonnors_XT& ndonors, Donnors_XT& donors,
                         Stack_XT& stack,
                         const Active_XT& active_nodes,
                         const Elevation_XT elevation,
                         typename Elevation_XT::value_type dx,
                         typename Elevation_XT::value_type dy)
{
    {
        PROFILE_COUNT(t0, "compute_basins", 8);
    basin_graph.compute_basins(basins, stack, receivers);
    //std::cout << basins << std::endl;
    }
    {
        PROFILE_COUNT(t1, "update_receivers", 8);

    basin_graph.template update_receivers<algo, connect>(receivers, dist2receivers, basins, stack,  active_nodes,
                                 elevation, dx, dy);
    }

    compute_donors(ndonors, donors, receivers);

    compute_stack(stack, ndonors, donors, receivers);
}

template <class Elevation_XT, class Stack_XT, class Rcv_XT>
void fill_sinks_flat(Elevation_XT& elevation, const Stack_XT& stack, const Rcv_XT& receivers)
{
    for(auto inode : stack)
    {
        auto ircv = receivers(inode);
        elevation(inode) = std::max(elevation(inode), elevation(ircv));
    }
}

template <class Water_XT, class Elevation_XT, class Stack_XT, class Rcv_XT>
void fill_sinks_flat(Water_XT& water, const Elevation_XT& elevation, const Stack_XT& stack, const Rcv_XT& receivers)
{
    for(auto inode : stack)
    {
        auto ircv = receivers(inode);
        if(inode == ircv)
            water(inode) = elevation(inode);
        else
            water(inode) = std::max(elevation(inode), water(ircv));
    }
}

template <class Elevation_XT, class Water_XT, class Active_XT>
bool check_fill_flat(const Elevation_XT& elevation, const Water_XT& water, const Active_XT& active_nodes)
{
    using elev_t = typename Elevation_XT::value_type;

    const auto elev_shape = elevation.shape();
    const index_t nrows = (index_t) elev_shape[0];
    const index_t ncols = (index_t) elev_shape[1];

    for(index_t r=0; r<nrows; ++r)
    {
        for(index_t c=0; c<ncols; ++c)
        {
            index_t inode = r * ncols + c;

            if(!active_nodes(inode))
                continue;

            if(elevation(inode) == water(inode))
                continue;

            for(size_t k=1; k<=8; ++k)
            {
                index_t kr = r + fs::consts::d8_row_offsets[k];
                index_t kc = c + fs::consts::d8_col_offsets[k];

                if(!fs::detail::in_bounds(elev_shape, kr, kc))
                    continue;

                index_t ineighbor = kr * ncols + kc;

                if(elevation(inode) > water(ineighbor))
                    return false;
            }
        }
    }
    return true;
}


template <BasinAlgo algo, ConnectType connect,
        class Elevation_XT, class Active_XT>
void fill_sinks_flat_basin_graph(Elevation_XT& elevation,
                                 const Active_XT& active_nodes,
                                 typename Elevation_XT::value_type dx,
                                 typename Elevation_XT::value_type dy)
{
    const auto elev_shape = elevation.shape();
    const size_t nrows = (size_t) elev_shape[0];
    const size_t ncols = (size_t) elev_shape[1];
    std::array<std::size_t, 1> shape_1D {nrows * ncols};

    BasinGraph<index_t, index_t, typename Elevation_XT::value_type> basin_graph;

    xt::xtensor<index_t, 1> basins (shape_1D);
    xt::xtensor<index_t, 1> receivers (shape_1D);
    xt::xtensor<typename Elevation_XT::value_type, 1> dist2receivers (shape_1D);
    xt::xtensor<index_t, 1> ndonors (shape_1D);
    xt::xtensor<index_t, 2> donors ({nrows * ncols, 8});
    xt::xtensor<index_t, 1> stack (shape_1D);

    {
        PROFILE(s0, "compute_receivers_d8");
        compute_receivers_d8(receivers, dist2receivers,
                             elevation, active_nodes,
                             1., 1.);
        //std::cout << receivers << std::endl;
    }


    {PROFILE(s1, "compute_donors");
        compute_donors(ndonors, donors, receivers);
    }

    {PROFILE(s2, "compute_stack");
        compute_stack(stack, ndonors, donors, receivers);

    }
    {
        PROFILE(s3, "correct_flowrouting");
        correct_flowrouting<algo, connect>(basin_graph, basins, receivers, dist2receivers,
                            ndonors, donors, stack,
                            active_nodes, elevation, dx, dy);

    }

    xt::xtensor<typename Elevation_XT::value_type, 1> water (shape_1D);

    fill_sinks_flat(water, elevation, stack, receivers);

    if(!check_fill_flat(elevation, water, active_nodes))
    {

        std::cout << elevation << std::endl;
        std::cout << water << std::endl;

        assert(false);
    }


    /*

    {
        PROFILE(s4, "fill_sinks_flat");
        fill_sinks_flat(elevation, stack, receivers);

    }*/
}



}  // namespace fastscapelib
