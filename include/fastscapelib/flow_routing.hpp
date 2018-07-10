/**
 * @brief Functions used to route (water) flow on a topographic
 * surface and compute flow path-related features or structures.
 *
 */
#pragma once

#include <cmath>
#include <algorithm>
#include <limits>
#include <array>
#include <type_traits>

#include "xtensor/xadapt.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xstrided_view.hpp"

#include "fastscapelib/utils.hpp"
#include "fastscapelib/consts.hpp"

#include "fastscapelib/basin_graph.hpp"

#include "fastscapelib/Profile.h"
#include <assert.h>
#include <xtensor/xio.hpp>

namespace fastscapelib
{

namespace detail
{


inline auto get_d8_distances(double dx, double dy) -> std::array<double, 9>
{
    std::array<double, 9> d8_dists;

    for(size_t k=0; k<9; ++k)
    {
        double d8_dx = dx * fastscapelib::consts::d8_col_offsets[k];
        double d8_dy = dy * fastscapelib::consts::d8_row_offsets[k];
        d8_dists[k] = std::sqrt(d8_dy*d8_dy + d8_dx*d8_dx);
    }

    return d8_dists;
}


inline auto get_d8_distances_inv(double dx, double dy) -> std::array<double, 9>
{
    std::array<double, 9> d8_dists;

    for(size_t k=0; k<9; ++k)
    {
        double d8_dx = dx * fastscapelib::consts::d8_col_offsets[k];
        double d8_dy = dy * fastscapelib::consts::d8_row_offsets[k];
        d8_dists[k] = 1.0/std::sqrt(d8_dy*d8_dy + d8_dx*d8_dx);
    }

    return d8_dists;
}


/**
 * compute_receivers_d8 implementation.
 */
template<class R, class D, class E, class A>
void compute_receivers_d8_impl(R&& receivers,
                               D&& dist2receivers,
                               E&& elevation,
                               A&& active_nodes,
                               double dx,
                               double dy)
{
    using elev_t = typename std::decay_t<E>::value_type;
    const auto d8_dists = detail::get_d8_distances(dx, dy);
    const auto d8_dists_inv = detail::get_d8_distances_inv(dx, dy);


    const auto elev_shape = elevation.shape();
    const auto nrows = static_cast<index_t>(elev_shape[0]);
    const auto ncols = static_cast<index_t>(elev_shape[1]);

    for(index_t r=0; r<nrows; ++r)
    {
        for(index_t c=0; c<ncols; ++c)
        {
            const index_t inode = r * ncols + c;

            receivers(inode) = inode;
            dist2receivers(inode) = 0.;

            if(!active_nodes(r, c))
            {
                continue;
            }

            double slope_max = 0.0;
            size_t k_found = 0;

            for(size_t k=1; k<=8; ++k)
            {
                const index_t kr = r + fastscapelib::consts::d8_row_offsets[k];
                const index_t kc = c + fastscapelib::consts::d8_col_offsets[k];

                if(!fastscapelib::detail::in_bounds(elev_shape, kr, kc))
                {
                    continue;
                }

                const index_t ineighbor = kr * ncols + kc;
                const double slope = (elevation(r, c) - elevation(kr, kc)) / d8_dists[k];

                if(slope > slope_max)
                {
                    slope_max = slope;
                    k_found = k;
                }
            }

            index_t kr = r + fastscapelib::consts::d8_row_offsets[k_found];
            index_t kc = c + fastscapelib::consts::d8_col_offsets[k_found];

            receivers(inode) = kr * ncols + kc;
            dist2receivers(inode) = d8_dists[k_found];
        }
    }
}


/**
 * compute_donors implementation.
 */
template<class N, class D, class R>
void compute_donors_impl(N&& ndonors,
                         D&& donors,
                         R&& receivers)
{
    const auto nnodes = static_cast<index_t>(receivers.size());

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


/**
 * non-recursive compute_stack implementation.
 */
template<class S, class N, class D, class R>
void compute_stack_impl(S&& stack,
                        N&& ndonors,
                        D&& donors,
                        R&& receivers)
{
    const auto nnodes = static_cast<index_t>(receivers.size());
    index_t nstack = 0;

    std::stack<index_t> tmp;

    for(index_t inode=0; inode<nnodes; ++inode)
    {
        if(receivers(inode) == inode)
        {
            tmp.push(inode);
            stack(nstack++) = inode;

            while(!tmp.empty())
            {
                index_t istack = tmp.top();
                tmp.pop();

                for(unsigned short k=0; k<ndonors(istack); ++k)
                {
                    index_t idonor = donors(istack, k);
                    stack(nstack++) = idonor;
                    tmp.push(idonor);
                }
            }
        }
    }
    assert(nstack == nnodes);
}


/**
 * compute_basins implementation.
 *
 * TODO: maybe return outlets as a new xtensor instead of nbasins?
 * TODO: maybe rename to `get_basin_ids` ? maybe not if return outlets
 */
template <class B, class O, class S, class R>
index_t compute_basins_impl(B&& basins,
                            O&& outlets_or_pits,
                            S&& stack,
                            R&& receivers)
{

    using Basin_T = typename std::decay_t<B>::value_type;
    using Node_T = typename std::decay_t<S>::value_type;

    fastscapelib::BasinGraph<Basin_T, Node_T> basin_graph;

    basin_graph.compute_basins(basins, stack, receivers);

    index_t nbasins = basin_graph.basin_count();
    auto bg_outlets = basin_graph.outlets();

    for (index_t i=0; i<nbasins; ++i)
    {
        outlets_or_pits(i) = bg_outlets[i];
    }

    return (nbasins);
}


/**
 * find_pits implementation.
 */
template<class P, class O, class A>
index_t find_pits_impl(P&& pits,
                       O&& outlets_or_pits,
                       A&& active_nodes,
                       index_t nbasins)
{
    index_t ipit = 0;
    const auto active_nodes_flat = xt::flatten(active_nodes);

    for(index_t ibasin=0; ibasin<nbasins; ++ibasin)
    {
        const index_t inode = outlets_or_pits(ibasin);

        if(active_nodes_flat(inode))
        {
            pits(ipit) = inode;
            ++ipit;
        }
    }

    index_t npits = ipit;

    return npits;
}


/**
 * compute_drainage_area implementation.
 */
template<class D, class C, class S, class R>
void compute_drainage_area_impl(D&& drainage_area,
                                C&& cell_area,
                                S&& stack,
                                R&& receivers)
{
    // reset drainage area values (must use a view to prevent resizing
    // drainage_area to 0-d when cell_area is 0-d!)
    auto drainage_area_ = xt::view(drainage_area, xt::all());
    drainage_area_ = cell_area;

    auto drainage_area_flat = xt::flatten(drainage_area_);

    for(auto inode=stack.crbegin(); inode!=stack.crend(); ++inode)
    {
        if(receivers(*inode) != *inode)
        {
            drainage_area_flat(receivers(*inode)) += drainage_area_flat(*inode);
        }
    }
}

}  // namespace detail


/**
 * Compute flow receivers on a rectangular grid using D8 single flow
 * routing method.
 *
 * Each node on the grid is assigned one receiver among its 8 direct
 * neighboring nodes according to the steepest slope (O’Callaghan and
 * Mark, 1984).
 *
 * When no downslope neighbor exist (i.e., pit nodes or grid
 * boundaries), the assigned receiver is the node itself. When two or
 * more neighbors have the same slope, the chosen neighbor is the
 * first one considered by the algorithm.
 *
 * This function also computes the planimetric distance between the
 * node and its receiver, which equals to the grid spacing in x or y
 * or to the distance between two diagonal neighbor nodes or 0.
 *
 * @param receivers : ``[intent=out, shape=(nnodes)]``
 *     Index of flow receiver at grid node.
 * @param dist2receivers : ``[intent=out, shape=(nnodes)]``
 *     Distance to receiver at grid node.
 * @param elevation : ``[intent=in, shape=(nrows, ncols)]``
 *     Topographic elevation at grid node
 * @param active_nodes : ``[intent=in, shape=(nrows, ncols)]``
 *     Boolean array for boundaries
 * @param dx : ``[intent=out]``
 *     Grid spacing in x
 * @param dy : ``[intent=out]``
 *     Grid spacing in y
 */
template<class R, class D, class E, class A>
void compute_receivers_d8(xtensor_t<R>& receivers,
                          xtensor_t<D>& dist2receivers,
                          const xtensor_t<E>& elevation,
                          const xtensor_t<A>& active_nodes,
                          double dx,
                          double dy)
{
    detail::compute_receivers_d8_impl(receivers.derived_cast(),
                                      dist2receivers.derived_cast(),
                                      elevation.derived_cast(),
                                      active_nodes.derived_cast(),
                                      dx, dy);
}


/**
 * Compute flow donors for each grid/mesh node.
 *
 * Flow donors are retrieved by simply inverting flow
 * receivers.
 *
 * @param ndonors : ``[intent=out, shape=(nnodes)]``
 *     Number of flow donors at grid node.
 * @param donors : ``[intent=out, shape=(nnodes, :)]``
 *     Indexes of flow donors at grid node.
 * @param receivers : ``[intent=in, shape=(nnodes)]``
 *     Index of flow receiver at grid node.
 */
template<class N, class D, class R>
void compute_donors(xtensor_t<N>& ndonors,
                    xtensor_t<D>& donors,
                    const xtensor_t<R>& receivers)
{
    detail::compute_donors_impl(ndonors.derived_cast(),
                                donors.derived_cast(),
                                receivers.derived_cast());
}


/**
 * Compute a stack of grid/mesh nodes to be used for flow tree
 * traversal.
 *
 * The stack is calculated recursively from outlets (or sinks) to
 * sources, using Braun and Willet's (2013) algorithm.
 *
 * @param stack : ``[intent=out, shape=(nnodes)]``
 *     Stack position at grid node.
 * @param ndonors : ``[intent=in, shape=(nnodes)]``
 *     Number of flow donors at grid node.
 * @param donors : ``[intent=in, shape=(nnodes, :)]``
 *     Indexes of flow donors at grid node.
 * @param receivers : ``[intent=in, shape=(nnodes)]``
 *     Index of flow receiver at grid node.
 */
template<class S, class N, class D, class R>
void compute_stack(xtensor_t<S>& stack,
                   const xtensor_t<N>& ndonors,
                   const xtensor_t<D>& donors,
                   const xtensor_t<R>& receivers)
{
    detail::compute_stack_impl(stack.derived_cast(),
                               ndonors.derived_cast(),
                               donors.derived_cast(),
                               receivers.derived_cast());
}


/**
 * Assign an id (integer) to each node of the grid/mesh that
 * corresponds to the catchment to which it belongs.
 *
 * A catchment (or drainage basin) is defined by an ensemble of
 * adjacent nodes through which all flow converges towards a common,
 * single node (outlet or pit).
 *
 * The algorithm performs a single traversal of the flow tree (in the
 * stack order) and increments the catchment id each time an outlet or
 * a pit is found (i.e., when the index of the flow receiver equals
 * the index of the node itself).
 *
 * This functions also computes the grid/mesh node indexes of
 * catchment outlets (or pits) and returns the total number of
 * catchments found inside the domain.
 *
 * @param basins: ``[intent=out, shape=(nnodes)]``
 *     Basin id at grid node.
 * @param outlets_or_pits : ``[intent=out, shape=(nnodes)]``
 *     Grid node index of the outlet (or pit)
 *     for basin id=0,1,...,nbasins-1.
 * @param stack :``[intent=in, shape=(nnodes)]``
 *     Stack position at grid node.
 * @param receivers : ``[intent=in, shape=(nnodes)]``
 *     Index of flow receiver at grid node.
 *
 * @returns
 *     Total number of drainage basins
 *     (``1 <= nbasins <= nnodes``).
 */
template<class B, class O, class S, class R>
index_t compute_basins(xtensor_t<B>& basins,
                       xtensor_t<O>& outlets_or_pits,
                       const xtensor_t<S>& stack,
                       const xtensor_t<R>& receivers)
{
    return detail::compute_basins_impl(basins.derived_cast(),
                                       outlets_or_pits.derived_cast(),
                                       stack.derived_cast(),
                                       receivers.derived_cast());
}


/**
 * Find grid/mesh nodes that are pits.
 *
 * @param pits: ``[intent=out, shape=(nnodes)]``
 *     Grid node index of the pit for pits in [0,1,...,npits-1].
 * @param outlets_or_pits : ``[intent=in, shape=(nnodes)]``
 *     Grid node index of the outlet (or pit)
 *     for basin id=0,1,...,nbasins-1.
 * @param active_nodes : ``[intent=in, shape=(nrows, ncols)]``
 *     Boolean array for boundaries
 * @param nbasins : ``[intent=in]``
 *     Total number of drainage basins (``1 <= nbasins <= nnodes``).
 *
 * @returns
 *     Total number of pits found (``0 <= npits <= nbasins``).
 */
template<class P, class O, class A>
index_t find_pits(xtensor_t<P>& pits,
                  const xtensor_t<O>& outlets_or_pits,
                  const xtensor_t<A>& active_nodes,
                  index_t nbasins)
{
    return detail::find_pits_impl(pits.derived_cast(),
                                  outlets_or_pits.derived_cast(),
                                  active_nodes.derived_cast(),
                                  nbasins);
}


/**
 * Compute drainage area for each node on a generic grid/mesh.
 *
 * @param drainage_area : ``[intent=out, shape=(nrows, ncols)||(nnodes)]``
 *     Drainage area at grid node.
 * @param cell_area : ``[intent=in, shape=(nrows, ncols)||(nnodes)||()]``
 *     Grid/mesh cell area at grid node (also accepts a scalar).
 * @param stack :``[intent=in, shape=(nnodes)]``
 *     Stack position at grid node.
 * @param receivers : ``[intent=in, shape=(nnodes)]``
 *     Index of flow receiver at grid node.
 */
template<class D, class C, class S, class R>
void compute_drainage_area(xtensor_t<D>& drainage_area,
                           const xtensor_t<C>& cell_area,
                           const xtensor_t<S>& stack,
                           const xtensor_t<R>& receivers)
{
    detail::compute_drainage_area_impl(drainage_area.derived_cast(),
                                       cell_area.derived_cast(),
                                       stack.derived_cast(),
                                       receivers.derived_cast());
}


/**
 * Compute drainage area for each node on a 2-d rectangular grid.
 *
 * @param drainage_area : ``[intent=inout, shape=(nrows, ncols)]``
 *     Drainage area at grid node.
 * @param stack :``[intent=in, shape=(nnodes)]``
 *     Stack position at grid node.
 * @param receivers : ``[intent=in, shape=(nnodes)]``
 *     Index of flow receiver at grid node.
 * @param dx : ``[intent=in]``
 *     Grid spacing in x.
 * @param dy : ``[intent=in]``
 *     Grid spacing in y.
 */
template<class D, class S, class R>
void compute_drainage_area(xtensor_t<D>& drainage_area,
                           const xtensor_t<S>& stack,
                           const xtensor_t<R>& receivers,
                           double dx,
                           double dy)
{
    xt::xtensor<double, 0> cell_area = dx * dy;
    compute_drainage_area(drainage_area, cell_area,  stack, receivers);
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
//    std::cout << elevation << std::endl;
//    std::cout << receivers << std::endl;

    {
        PROFILE_COUNT(t0, "compute_basins", 8);
    basin_graph.compute_basins(basins, stack, receivers);
//    std::cout << basins << std::endl;
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

template <class Elevation_XT, class Water_XT, class Active_XT, class Rcv_XT>
bool check_fill_flat(const Elevation_XT& elevation, const Water_XT& water, const Active_XT& active_nodes, const Rcv_XT& rcv)
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

            if(rcv(inode) == inode)
            {
                std::cout << "Pits remaining !" << std::endl;
                return false;
            }

            if(elevation(inode) == water(inode))
                continue;

            for(size_t k=1; k<=8; ++k)
            {
                index_t kr = r + fastscapelib::consts::d8_row_offsets[k];
                index_t kc = c + fastscapelib::consts::d8_col_offsets[k];

                if(!fastscapelib::detail::in_bounds(elev_shape, kr, kc))
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
                             dx, dy);
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

    if(!check_fill_flat(elevation, water, active_nodes, receivers))
    {

        std::cout << basins << std::endl;
        std::cout << elevation << std::endl;
        std::cout << water << std::endl;

        assert(false);
        exit(0);
    }


    /*

    {
        PROFILE(s4, "fill_sinks_flat");
        fill_sinks_flat(elevation, stack, receivers);

    }*/
}



}  // namespace fastscapelib
