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


namespace fs = fastscapelib;


namespace fastscapelib
{

namespace detail
{

    inline auto get_d8_distances(double dx, double dy) -> std::array<double, 9>
    {
        std::array<double, 9> d8_dists;

        for(size_t k=0; k<9; ++k)
        {
            d8_dists[k] = std::sqrt(
                std::pow(dy * fs::consts::d8_row_offsets[k], 2.0) +
                std::pow(dx * fs::consts::d8_col_offsets[k], 2.0));
        }

        return d8_dists;
    }


    template<class D1, class D2, class D3>
    void add2stack(index_t& nstack,
                   D1& stack,
                   const D2& ndonors,
                   const D3& donors,
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
    R& d_receivers = receivers.derived_cast();
    D& d_dist2receivers = dist2receivers.derived_cast();
    const E& d_elevation = elevation.derived_cast();
    const A& d_active_nodes = active_nodes.derived_cast();

    using elev_t = typename E::value_type;
    const auto d8_dists = detail::get_d8_distances(dx, dy);

    const auto elev_shape = d_elevation.shape();
    const index_t nrows = (index_t) elev_shape[0];
    const index_t ncols = (index_t) elev_shape[1];

    for(index_t r=0; r<nrows; ++r)
    {
        for(index_t c=0; c<ncols; ++c)
        {
            index_t inode = r * ncols + c;

            d_receivers(inode) = inode;
            d_dist2receivers(inode) = 0.;

            if(!d_active_nodes(r, c))
            {
                continue;
            }

            double slope_max = std::numeric_limits<double>::min();

            for(size_t k=1; k<=8; ++k)
            {
                index_t kr = r + fs::consts::d8_row_offsets[k];
                index_t kc = c + fs::consts::d8_col_offsets[k];

                if(!fs::detail::in_bounds(elev_shape, kr, kc))
                {
                    continue;
                }

                index_t ineighbor = kr * ncols + kc;
                double slope = (d_elevation(r, c) - d_elevation(kr, kc)) / d8_dists[k];

                if(slope > slope_max)
                {
                    slope_max = slope;
                    d_receivers(inode) = ineighbor;
                    d_dist2receivers(inode) = d8_dists[k];
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
 * @param[out]  ndonors    Number of flow donors at grid node.
 *                           ``[shape=(nnodes)]``
 * @param[out]  donors     Indexes of flow donors at grid node.
 *                           ``[shape=(nnodes, :)]``
 * @param[in]   receivers  Index of flow receiver at grid node.
 *                           ``[shape=(nnodes)]``
 */
template<class Xndonors, class Xdonors, class Xrec>
void compute_donors(Xndonors& ndonors,
                    Xdonors& donors,
                    const Xrec& receivers)
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
    S& d_stack = stack.derived_cast();
    const N& d_ndonors = ndonors.derived_cast();
    const D& d_donors = donors.derived_cast();
    const R& d_receivers = receivers.derived_cast();

    index_t nnodes = (index_t) d_receivers.size();
    index_t nstack = 0;

    for(index_t inode=0; inode<nnodes; ++inode)
    {
        if(d_receivers(inode) == inode)
        {
            d_stack(nstack) = inode;
            ++nstack;
            detail::add2stack(nstack, d_stack, d_ndonors, d_donors, inode);
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
 * @param[out]  basins           Basin id at grid node.
 *                                 ``[shape=(nnodes)]``
 * @param[out]  outlets_or_pits  Grid node index of the outlet (or pit)
 *                               for basin id=0,1,...,nbasins-1.
 *                                 ``[shape=(nnodes)]``
 * @param[out]  stack            Stack position at grid node.
 *                                 ``[shape=(nnodes)]``
 * @param[in]   receivers        Index of flow receiver at grid node.
 *                                 ``[shape=(nnodes)]``
 *
 * @returns  Total number of drainage basins
 *           (``1 <= nbasins <= nnodes``).
 */
template<class Xbasins, class Xoutlets, class Xstack, class Xrec>
index_t compute_basins(Xbasins& basins,
                       Xoutlets& outlets_or_pits,
                       const Xstack& stack,
                       const Xrec& receivers)
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


template<class A1, class A2, class A3>
index_t compute_pits(A1& pits,
                     const A2& outlets,
                     const A3& active_nodes,
                     index_t nbasins)
{
    //TODO: works whether active_nodes is 1-d or 2-d but not safe!
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

}  // namespace fastscapelib
