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
#include "xtensor/xstrided_view.hpp"

#include "fastscapelib/utils.hpp"
#include "fastscapelib/consts.hpp"


namespace fs = fastscapelib;


namespace fastscapelib
{

namespace detail
{

    inline auto get_d8_distances(double dx, double dy) -> std::array<double, 9>
    {
        std::array<double, 9> d8_dists;

        for(unsigned short k=0; k<9; ++k)
        {
            d8_dists[k] = std::sqrt(
                std::pow(dy * fs::consts::d8_row_offsets[k], 2.0) +
                std::pow(dx * fs::consts::d8_col_offsets[k], 2.0));
        }

        return d8_dists;
    }


    template<class S, class N, class D>
    void add2stack(index_t& nstack,
                   xtensor_t<S>& stack,
                   const xtensor_t<N>& ndonors,
                   const xtensor_t<D>& donors,
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
    using elev_t = typename E::value_type;
    const auto d8_dists = detail::get_d8_distances(dx, dy);

    const auto elev_shape = elevation.shape();
    const index_t nrows = static_cast<index_t>(elev_shape[0]);
    const index_t ncols = static_cast<index_t>(elev_shape[1]);

    for(index_t r=0; r<nrows; ++r)
    {
        for(index_t c=0; c<ncols; ++c)
        {
            index_t inode = r * ncols + c;

            receivers(inode) = inode;
            dist2receivers(inode) = 0.;

            if(!active_nodes(r, c))
            {
                continue;
            }

            double slope_max = std::numeric_limits<double>::min();

            for(unsigned short k=1; k<=8; ++k)
            {
                index_t kr = r + fs::consts::d8_row_offsets[k];
                index_t kc = c + fs::consts::d8_col_offsets[k];

                if(!fs::detail::in_bounds(elev_shape, kr, kc))
                {
                    continue;
                }

                index_t ineighbor = kr * ncols + kc;
                double slope = (elevation(r, c) - elevation(kr, kc)) / d8_dists[k];

                if(slope > slope_max)
                {
                    slope_max = slope;
                    receivers(inode) = ineighbor;
                    dist2receivers(inode) = d8_dists[k];
                }
            }
        }
    }
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
    index_t nnodes = static_cast<index_t>(receivers.size());

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
    index_t nnodes = static_cast<index_t>(receivers.size());
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
    index_t ibasin = -1;

    for(auto&& istack : stack)
    {
        index_t irec = receivers(istack);

        if(irec == istack)
        {
            ++ibasin;
            outlets_or_pits(ibasin) = istack;
        }

        basins(istack) = ibasin;
    }

    index_t nbasins = ibasin + 1;

    return nbasins;
}


template<class P, class O, class A>
index_t compute_pits(xtensor_t<P>& pits,
                     const xtensor_t<O>& outlets,
                     const xt::xexpression<A>& active_nodes,
                     index_t nbasins)
{
    index_t ipit = 0;
    const A& active_nodes_ = active_nodes.derived_cast();
    const auto active_nodes_flat = xt::flatten(active_nodes_);

    for(index_t ibasin=0; ibasin<nbasins; ++ibasin)
    {
        index_t inode = outlets(ibasin);

        if(active_nodes_flat(inode))
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

    auto area_flat = xt::flatten(area);

    compute_drainage_area(area_flat, stack, receivers);
}

}  // namespace fastscapelib
