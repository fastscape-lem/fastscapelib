/**
 * Grid boundaries definitions and helper functions.
 */
#pragma once

#include <cstdint>
#include <map>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xview.hpp"

#include "fastscapelib/utils.hpp"


namespace fastscapelib
{


enum class NodeStatus : std::int8_t
{
    CORE_NODE = 0,
    FIXED_VALUE_BOUNDARY = 1,
    FIXED_GRADIENT_BOUNDARY = 2,
    LOOPED_BOUNDARY = 3
};


template<class S>
auto create_node_status(S& shape)
{
    return xt::zeros<NodeStatus>(shape);
}


template<class NS>
void set_node_status_grid_boundaries(NS& node_status,
                                     NodeStatus top,
                                     NodeStatus right,
                                     NodeStatus bottom,
                                     NodeStatus left)
{
    const auto shape = node_status.shape();
    auto last_row = static_cast<index_t>(shape[0] - 1);
    auto last_col = static_cast<index_t>(shape[1] - 1);

    auto top_nodes = xt::view(node_status, 0, xt::all());
    top_nodes = top;

    auto right_nodes = xt::view(node_status, xt::all(), last_col);
    right_nodes = right;

    auto bottom_nodes = xt::view(node_status, last_row, xt::all());
    bottom_nodes = bottom;

    auto left_nodes = xt::view(node_status, xt::all(), 0);
    left_nodes = left;

    // maybe fix grid corners
    std::map<NodeStatus, std::int8_t> priority = {
        {NodeStatus::CORE_NODE, 0},
        {NodeStatus::LOOPED_BOUNDARY, 1},
        {NodeStatus::FIXED_GRADIENT_BOUNDARY, 2},
        {NodeStatus::FIXED_VALUE_BOUNDARY, 3}
    };

    if (priority[right]  < priority[top])    node_status(0, last_col) = top;
    if (priority[bottom] < priority[right])  node_status(last_row, last_col) = right;
    if (priority[left]   < priority[bottom]) node_status(last_row, 0) = bottom;
    if (priority[left]   < priority[top])    node_status(0, 0) = top;

}

}  // namespace fastscapelib
