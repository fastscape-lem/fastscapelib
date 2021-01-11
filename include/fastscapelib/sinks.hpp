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
        template<class T>
        struct node_container
        {
            index_t r;
            index_t c;
            T value;
            node_container() { }
            node_container(index_t row, index_t col, T val) :
                r(row), c(col), value(val) { }
            bool operator > (const node_container<T>& other) const
            {
                return value > other.value;
            }
        };


        template<class T>
        using node_pr_queue = std::priority_queue<node_container<T>,
                                                std::vector<node_container<T>>,
                                                std::greater<node_container<T>>>;


        template<class T>
        using node_queue = std::queue<node_container<T>>;


        /**
         * Initialize priority flood algorithms.
         *
         * Add border grid nodes to the priority queue and mark them as
         * resolved.
         */
        template<class E, class elev_t = typename std::decay_t<E>::value_type>
        void init_pflood(E&& elevation,
                        xt::xtensor<bool, 2>& closed,
                        node_pr_queue<elev_t>& open)
        {
            auto elev_shape = elevation.shape();

            index_t nrows = static_cast<index_t>(elev_shape[0]);
            index_t ncols = static_cast<index_t>(elev_shape[1]);

            auto place_node = [&](index_t row, index_t col)
            {
                open.emplace(node_container<elev_t>(row, col, elevation(row, col)));
                closed(row, col) = true;
            };

            for(index_t c=0; c<ncols; ++c)
            {
                place_node(0, c);
                place_node(nrows-1, c);
            }

            for(index_t r=1; r<nrows-1; ++r)
            {
                place_node(r, 0);
                place_node(r, ncols-1);
            }
        }


        /**
         * fill_sinks_flat implementation.
         */
        template<class E>
        void fill_sinks_flat_impl(E&& elevation)
        {
            using elev_t = typename std::decay_t<E>::value_type;
            auto elev_shape = elevation.shape();

            node_pr_queue<elev_t> open;
            xt::xtensor<bool, 2> closed = xt::zeros<bool>(elev_shape);

            init_pflood(elevation, closed, open);

            while(open.size()>0)
            {
                node_container<elev_t> inode = open.top();
                open.pop();

                for(unsigned short k=1; k<=8; ++k)
                {
                    index_t kr = inode.r + fastscapelib::consts::d8_row_offsets[k];
                    index_t kc = inode.c + fastscapelib::consts::d8_col_offsets[k];

                    if(!fastscapelib::detail::in_bounds(elev_shape, kr, kc)) { continue; }
                    if(closed(kr, kc)) { continue; }

                    elevation(kr, kc) = std::max(elevation(kr, kc), inode.value);
                    open.emplace(node_container<elev_t>(kr, kc, elevation(kr, kc)));
                    closed(kr, kc) = true;
                }
            }
        }


        /**
         * fill_sinks_sloped implementation.
         */
        template<class E>
        void fill_sinks_sloped_impl(E&& elevation)
        {
            using elev_t = typename std::decay_t<E>::value_type;
            auto elev_shape = elevation.shape();

            node_pr_queue<elev_t> open;
            node_queue<elev_t> pit;
            xt::xtensor<bool, 2> closed = xt::zeros<bool>(elev_shape);

            init_pflood(elevation, closed, open);

            while(!open.empty() || !pit.empty())
            {
                node_container<elev_t> inode, knode;

                if(!pit.empty() &&
                (open.empty() || open.top().value == pit.front().value))
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
                    inode.value, std::numeric_limits<elev_t>::infinity());

                for(unsigned short k=1; k<=8; ++k)
                {
                    index_t kr = inode.r + fastscapelib::consts::d8_row_offsets[k];
                    index_t kc = inode.c + fastscapelib::consts::d8_col_offsets[k];

                    if(!fastscapelib::detail::in_bounds(elev_shape, kr, kc)) { continue; }
                    if(closed(kr, kc)) { continue; }

                    if(elevation(kr, kc) <= elev_tiny_step)
                    {
                        elevation(kr, kc) = elev_tiny_step;
                        knode = node_container<elev_t>(kr, kc, elevation(kr, kc));
                        pit.emplace(knode);
                    }
                    else
                    {
                        knode = node_container<elev_t>(kr, kc, elevation(kr, kc));
                        open.emplace(knode);
                    }

                    closed(kr, kc) = true;
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
    template<class E>
    void fill_sinks_flat(xtensor_t<E>& elevation)
    {
        detail::fill_sinks_flat_impl(elevation.derived_cast());
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
    template<class E>
    void fill_sinks_sloped(xtensor_t<E>& elevation)
    {
        detail::fill_sinks_sloped_impl(elevation.derived_cast());
    }
}  // namespace fastscapelib

#endif