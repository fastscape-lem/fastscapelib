/**
 * Grid boundaries definitions and helper functions.
 */
#pragma once

#include <cstdint>
#include <map>
#include <stdexcept>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xview.hpp"

#include "fastscapelib/utils.hpp"


namespace fastscapelib
{


/**
 * Status of grid/mesh nodes either inside the domain or on the domain
 * boundary.
 */
enum class NodeStatus : std::int8_t
{
    CORE_NODE = 0,
    FIXED_VALUE_BOUNDARY = 1,
    FIXED_GRADIENT_BOUNDARY = 2,
    LOOPED_BOUNDARY = 3
};


/**
 * Create an array and set all of its elements to
 * ``NodeStatus::CORE_NODE``.
 *
 * This helper function ensures that the array values for node status
 * have the right type.
 *
 * @param shape : ``[intent=in]``
 *     The shape of the array to create.
 *
 * @returns
 *     A new, uniform array with ``NodeStatus::CORE_NODE`` values.
 */
template<class S>
auto create_node_status(S& shape)
{
    return xt::zeros<NodeStatus>(shape);
}


/**
 * Helper function for setting node status at each side of the grid
 * boundaries.
 *
 * Conflicts may arise on the grid corners if different status are set
 * on adjacent grid sides. Those conflicts are solved using the
 * following priority order:
 *
 * FIXED_VALUE_BOUNDARY > FIXED_GRADIENT_BOUNDARY > LOOPED_BOUNDARY > CORE_NODE
 *
 * @param node_status : ``[intent=inout, shape=(nrows, ncols)]``
 *     Array of node status on the grid.
 * @param top :
 *     Node status to set on the top boundary (first row).
 * @param right :
 *     Node status to set on the right boundary (last col).
 * @param bottom :
 *     Node status to set on the bottom boundary (last row).
 * @param left :
 *     Node status to set on the left boundary (first col).
 *
 * @throws std::invalid_argument
 *     If looped boundary conditions are not symmetrical.
 */
template<class NS>
void set_node_status_grid_boundaries(NS& node_status,
                                     NodeStatus top,
                                     NodeStatus right,
                                     NodeStatus bottom,
                                     NodeStatus left)
{
    // check symmetry of looped boundaries
    if ((top == NodeStatus::LOOPED_BOUNDARY ^ bottom == NodeStatus::LOOPED_BOUNDARY) ||
        (left == NodeStatus::LOOPED_BOUNDARY ^ right == NodeStatus::LOOPED_BOUNDARY))
    {
        throw std::invalid_argument("looped boundaries are not symmetrical");
    }

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
