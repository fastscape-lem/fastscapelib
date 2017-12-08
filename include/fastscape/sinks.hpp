/**
 * @file
 * @brief Provides implementation of various efficient
 *        algorithms for depression filling or pit resolving.
 */

#pragma once

#include <cmath>
#include <functional>
#include <algorithm>
#include <queue>
#include <limits>
#include <stdexcept>


#include "xtensor/xtensor.hpp"

#include "fastscape/utils.hpp"
#include "fastscape/consts.hpp"


namespace fs = fastscape;


namespace fastscape {

    namespace detail {

        /**
         * @brief A simple grid node container.
         *
         * Stores both a position (r, c) and a value at that position.
         * Also defines  operator '>' that compares only on `value`.
         *
         * The main purpose of this container is for using with
         * (priority) queues.
         */
        template<class T>
        struct node_container {
            index_t r;
            index_t c;
            T value;
            node_container() { }
            node_container(index_t row, index_t col, T val) :
                r(row), c(col), value(val) { }
            bool operator > (const node_container<T>& other) const {
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
         * @brief Return true if a given (row, col) index is in bounds (false
         *        otherwise).
         */
        template<class S>
        bool in_bounds(const S& shape, index_t row, index_t col) {
            if(row >= 0 && row < (index_t) shape[0]
                   && col >= 0 && col < (index_t) shape[1]) {
                return true;
            }
            return false;
        }

        /**
         * @brief Initialize priority flood algorithms.
         *
         * Add border grid nodes to the priority queue and mark them as
         * resolved.
         */
        template<class E, class elev_t>
        void init_pflood(xt::xexpression<E>& elevation,
                         xt::xtensor<bool, 2>& closed,
                         node_pr_queue<elev_t>& open) {
            auto& elev = elevation.derived_cast();
            auto elev_shape = elev.shape();

            index_t nrows = (index_t) elev_shape[0];
            index_t ncols = (index_t) elev_shape[1];

            auto place_node = [&](index_t row, index_t col) {
                open.emplace(node_container<elev_t>(row, col, elev(row, col)));
                closed(row, col) = true;
            };

            for(index_t c=0; c<ncols; ++c) {
                place_node(0, c);
                place_node(nrows-1, c);
            }

            for(index_t r=1; r<nrows-1; ++r) {
                place_node(r, 0);
                place_node(r, ncols-1);
            }
        }

    }

    template<class E>
    auto fill_sinks_flat(xt::xexpression<E>& elevation) {
        using elev_t = typename E::value_type;
        auto& elev = elevation.derived_cast();
        auto elev_shape = elev.shape();

        detail::node_pr_queue<elev_t> open;
        xt::xtensor<bool, 2> closed = xt::zeros<bool>(elev_shape);

        detail::init_pflood(elevation, closed, open);

        while(open.size()>0) {
            detail::node_container<elev_t> inode = open.top();
            open.pop();

            for(int k=1; k<=8; ++k) {
                index_t kr = inode.r + fs::consts::d8_row_offsets[k];
                index_t kc = inode.c + fs::consts::d8_col_offsets[k];

                if(!detail::in_bounds(elev_shape, kr, kc)) { continue; }
                if(closed(kr, kc)) { continue; }

                elev(kr, kc) = std::max(elev(kr, kc), inode.value);
                open.emplace(detail::node_container<elev_t>(kr, kc, elev(kr, kc)));
                closed(kr, kc) = true;
            }
        }
    }


    template<class E>
        auto fill_sinks_sloped(xt::xexpression<E>& elevation) {
        using elev_t = typename E::value_type;
        auto& elev = elevation.derived_cast();
        auto elev_shape = elev.shape();

        detail::node_pr_queue<elev_t> open;
        detail::node_queue<elev_t> pit;
        xt::xtensor<bool, 2> closed = xt::zeros<bool>(elev_shape);

        detail::init_pflood(elevation, closed, open);

        while(!open.empty() || !pit.empty()) {
            detail::node_container<elev_t> inode, knode;

            if(!pit.empty() &&
                   (open.empty() || open.top().value == pit.front().value)) {
                inode = pit.front();
                pit.pop();
            }
            else {
                inode = open.top();
                open.pop();
            }

            elev_t elev_tiny_step = std::nextafter(
                inode.value, std::numeric_limits<elev_t>::infinity());

            for(int k=1; k<=8; ++k) {
                index_t kr = inode.r + fs::consts::d8_row_offsets[k];
                index_t kc = inode.c + fs::consts::d8_col_offsets[k];

                if(!detail::in_bounds(elev_shape, kr, kc)) { continue; }
                if(closed(kr, kc)) { continue; }

                if(elev(kr, kc) <= elev_tiny_step) {
                    elev(kr, kc) = elev_tiny_step;
                    knode = detail::node_container<elev_t>(kr, kc, elev(kr, kc));
                    pit.emplace(knode);
                }
                else {
                    knode = detail::node_container<elev_t>(kr, kc, elev(kr, kc));
                    open.emplace(knode);
                }

                closed(kr, kc) = true;
            }
        }
    }

}
