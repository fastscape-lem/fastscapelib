/**
 * Provides implementation of various efficient algorithms for
 * depression filling or pit resolving.
 */
#ifndef FASTSCAPELIB_SINKS_H
#define FASTSCAPELIB_SINKS_H

#include <cmath>
#include <functional>
#include <algorithm>
#include <queue>
#include <limits>
#include <type_traits>

#include "xtensor/xtensor.hpp"

#include "fastscapelib/utils.hpp"
#include "fastscapelib/consts.hpp"
#include "fastscapelib/structured_grid.hpp"


namespace fastscapelib
{
    namespace detail
    {


        /**
         * A simple grid node container.
         *
         * Stores both a position (r, c) and a value at that position.
         * Also defines  operator '>' that compares only on `value`.
         *
         * The main purpose of this container is for using with
         * (priority) queues.
         */
        template<class G, class T>
        struct node_container
        {
            using size_type = typename G::size_type;

            size_type m_idx;
            T m_elevation;

            node_container() { }
            node_container(size_type idx, T elevation) :
                m_idx(idx), m_elevation(elevation) { }
            bool operator > (const node_container<G, T>& other) const
            {
                return m_elevation > other.m_elevation;
            }
        };


        template<class G, class T>
        using node_pr_queue = std::priority_queue<node_container<G, T>,
                                                std::vector<node_container<G, T>>,
                                                std::greater<node_container<G, T>>>;


        template<class G, class T>
        using node_queue = std::queue<node_container<G, T>>;


        /**
         * Initialize priority flood algorithms.
         *
         * Add fixed value grid nodes to the priority queue and mark them as
         * resolved.
         */
        template<class G, class E, class elev_t = typename std::decay_t<E>::value_type>
        void init_pflood(G& grid,
                         E&& elevation,
                         xt::xtensor<bool, 1>& closed,
                         node_pr_queue<G, elev_t>& open)
        {
            using size_type = typename G::size_type;

            // TODO: assert elevation shape match grid shape

            const auto elevation_flat = xt::flatten(elevation);

            for (size_type idx=0; idx<grid.size(); ++idx)
            {
                if (grid.status_at_nodes()[idx] == fastscapelib::node_status::fixed_value_boundary)
                {
                    open.emplace(node_container<G, elev_t>(idx, elevation_flat(idx)));
                    closed(idx) = true;
                }
            }
        }


        /**
         * fill_sinks_flat implementation.
         */
        template<class G, class E>
        void fill_sinks_flat_impl(G& grid, E&& elevation)
        {
            using size_type = typename G::size_type;
            using elev_t = typename std::decay_t<E>::value_type;

            auto elevation_flat = xt::flatten(elevation);

            node_pr_queue<G, elev_t> open;
            xt::xtensor<bool, 1> closed = xt::zeros<bool>({grid.size()});

            init_pflood(grid, elevation, closed, open);

            while(open.size()>0)
            {
                node_container<G, elev_t> inode = open.top();
                open.pop();

                for(auto nidx : grid.neighbors_indices(inode.m_idx))
                {
                    if(closed(nidx))
                    {
                        continue;
                    }


                    elevation_flat(nidx) = std::max(elevation_flat(nidx), inode.m_elevation);
                    open.emplace(node_container<G, elev_t>(nidx, elevation_flat(nidx)));
                    closed(nidx) = true;
                }
            }
        }


        /**
         * fill_sinks_sloped implementation.
         */
        template<class G, class E>
        void fill_sinks_sloped_impl(G& grid, E&& elevation)
        {
            using size_type = typename G::size_type;
            using elev_t = typename std::decay_t<E>::value_type;

            auto elevation_flat = xt::flatten(elevation);

            node_pr_queue<G, elev_t> open;
            node_queue<G, elev_t> pit;
            xt::xtensor<bool, 1> closed = xt::zeros<bool>({grid.size()});

            init_pflood(grid, elevation, closed, open);

            while(!open.empty() || !pit.empty())
            {
                node_container<G, elev_t> inode, knode;

                if(!pit.empty() && (open.empty() || open.top().m_elevation == pit.front().m_elevation))
                {
                    inode = pit.front();
                    pit.pop();
                }
                else
                {
                    inode = open.top();
                    open.pop();
                }

                elev_t elev_tiny_step = std::nextafter(
                    inode.m_elevation, std::numeric_limits<elev_t>::infinity());

                for(auto&& nidx : grid.neighbors_indices(inode.m_idx))
                {
                    if(closed(nidx))
                    {
                        continue;
                    }

                    if(elevation_flat(nidx) <= elev_tiny_step)
                    {
                        elevation_flat(nidx) = elev_tiny_step;
                        knode = node_container<G, elev_t>(nidx, elevation_flat(nidx));
                        pit.emplace(knode);
                    }
                    else
                    {
                        knode = node_container<G, elev_t>(nidx, elevation_flat(nidx));
                        open.emplace(knode);
                    }

                    closed(nidx) = true;
                }
            }
        }

    }  // namespace detail


    /**
     * Fill all pits and remove all digital dams from a digital
     * elevation model (rectangular grid).
     *
     * Elevation values may be updated so that no depression (sinks or
     * local minima) remains.
     *
     * The algorithm is based on priority queues and is detailed
     * in Barnes (2014). It fills sinks with flat areas.
     *
     * @param elevation : ``[intent=inout, shape=(nrows, ncols)]``
     *     Elevation at grid node.
     */
    template<class G, class E>
    void fill_sinks_flat(G& grid, xtensor_t<E>& elevation)
    {
        detail::fill_sinks_flat_impl(grid, elevation.derived_cast());
    }


    /**
     * Fill all pits and remove all digital dams from a digital
     * elevation model (rectangular grid).
     *
     * Elevation values may be updated so that no depression (sinks or
     * local minima) remains.
     *
     * The algorithm is based on priority queues and is detailed in Barnes
     * (2014). This variant fills sinks with nearly flat areas
     * (i.e. elevation is increased by small amount) so that there is no
     * drainage singularities.
     *
     * @param elevation : ``[intent=inout, shape=(nrows, ncols)]``
     *     Elevation at grid node.
     *
     * @sa fill_sinks_flat
     */
    template<class G, class E>
    void fill_sinks_sloped(G& grid, xtensor_t<E>& elevation)
    {
        detail::fill_sinks_sloped_impl(grid, elevation.derived_cast());
    }
}  // namespace fastscapelib

#endif
