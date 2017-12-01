/**
 * @file
 * @brief Provides implementation of various efficient
 *        algorithms for depression filling.
 */

#ifndef DEPRESSIONS_HPP
#define DEPRESSIONS_HPP

#include <cmath>
#include <functional>
#include <algorithm>
#include <queue>

#include "xtensor/xtensor.hpp"

#include "constants.hpp"


/**
    @brief A simple grid node container.

    Stores both a position (r, c) and a value at that position.
    Also defines  operator '>' that compares only on `value`.

    The main purpose of this container is for using with (priority) queues.
 */
template<class T>
struct node_container {
    ssize_t r;
    ssize_t c;
    T value;
    node_container() { }
    node_container(ssize_t row, ssize_t col, T val) : r(row), c(col), value(val) { }
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


template<class E>
auto priority_flood_original(xt::xexpression<E>& elevation) {
    using elev_t = typename E::value_type;
    auto& elev = elevation.derived_cast();
    auto elev_shape = elev.shape();

    ssize_t nrows = elev_shape[0];
    ssize_t ncols = elev_shape[1];

    node_pr_queue<elev_t> open;
    xt::xtensor<bool, 2> closed = xt::zeros<bool>(elev_shape);

    auto place_node = [&](ssize_t r, ssize_t c) {
        open.emplace(node_container<elev_t>(r, c, elev(r, c)));
        closed(r, c) = true;
    };

    // Add top and bottom rows in priority queue
    for(ssize_t c=0; c<ncols; c++) {
        place_node(0, c);
        place_node(nrows-1, c);
    }

    // Add left and right borders in priority queue
    for(ssize_t r=1; r<nrows-1; r++) {
        place_node(r, 0);
        place_node(r, ncols-1);
    }

    while(open.size()>0) {
        node_container<elev_t> inode = open.top();
        open.pop();

        for(int ineighbor=1; ineighbor<=8; ineighbor++) {
            ssize_t inb_r = inode.r + fs::constants::d8_row_offsets[ineighbor];
            ssize_t inb_c = inode.c + fs::constants::d8_col_offsets[ineighbor];

            try { if(closed.at(inb_r, inb_c)) continue; }
            catch (const std::out_of_range& e) { continue; }

            elev(inb_r, inb_c) = std::max(elev(inb_r, inb_c), inode.value);
            place_node(inb_r, inb_c);
        }
    }
}


#endif
