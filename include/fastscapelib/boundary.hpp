/**
 * Grid boundaries definitions and helper functions.
 */
#pragma once

#include <cstdint>
#include <map>
#include <stdexcept>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xview.hpp"

#include "fastscapelib/grid.hpp"
#include "fastscapelib/xtensor_utils.hpp"
#include "fastscapelib/utils.hpp"


namespace fastscapelib
{


/**
 * Status of grid/mesh nodes either inside the domain or on the domain
 * boundary.
 */
enum class NodeStatus : std::uint8_t
{
    CORE_NODE = 0,
    FIXED_VALUE_BOUNDARY = 1,
    FIXED_GRADIENT_BOUNDARY = 2,
    LOOPED_BOUNDARY = 3
};


namespace detail
{


template <grid_type G, array_type A>
class grid_boundary_base
{
protected:

    using T = std::underlying_type_t<NodeStatus>;
    using X = xt_container_t<A, T, grid_properties<G>::ndim_nodes>;

public:

    X node_status;

    template <class S>
    grid_boundary_base(S& shape)
    {
        node_status = xt::zeros<T>(shape);
    }

    grid_boundary_base(std::initializer_list<std::size_t> shape)
    {
        node_status = xt::zeros<T>(shape);
    }

    ~grid_boundary_base() = default;
};


} // namespace detail


template <grid_type G, array_type A = array_type::xtensor>
class grid_boundary : public detail::grid_boundary_base<G, A>
{
public:

    template <class S>
    grid_boundary(S& shape) : detail::grid_boundary_base<G, A>(shape)
    {
    }

    grid_boundary(std::initializer_list<std::size_t> shape)
        : detail::grid_boundary_base<G, A>(shape)
    {
    }
};


template <array_type A>
class grid_boundary<grid_type::raster, A>
    : public detail::grid_boundary_base<grid_type::raster, A>
{
public:

    NodeStatus top = NodeStatus::CORE_NODE;
    NodeStatus right = NodeStatus::CORE_NODE;
    NodeStatus bottom = NodeStatus::CORE_NODE;
    NodeStatus left = NodeStatus::CORE_NODE;

    template <class S>
    grid_boundary(S& shape)
        : detail::grid_boundary_base<grid_type::raster, A>(shape)
    {
    }

    grid_boundary(std::initializer_list<std::size_t> shape)
        : detail::grid_boundary_base<grid_type::raster, A>(shape)
    {
    }

    void set_status_at_boundary(const NodeStatus border);

    void set_status_at_boundary(const NodeStatus top,
                                const NodeStatus right,
                                const NodeStatus bottom,
                                const NodeStatus left);

private:

    using T = typename detail::grid_boundary_base<grid_type::raster, A>::T;

    std::map<NodeStatus, T> priority
    {
        {NodeStatus::CORE_NODE, 0},
        {NodeStatus::LOOPED_BOUNDARY, 1},
        {NodeStatus::FIXED_GRADIENT_BOUNDARY, 2},
        {NodeStatus::FIXED_VALUE_BOUNDARY, 3}
    };
};


template <array_type A>
void grid_boundary<grid_type::raster, A>::set_status_at_boundary(const NodeStatus border)
{
    auto row_bound = xt::view(this->node_status, xt::keep(0, -1), xt::all());
    row_bound = utils::cast_underlying(border);

    auto col_bound = xt::view(this->node_status, xt::all(), xt::keep(0, -1));
    col_bound = utils::cast_underlying(border);
}


template <array_type A>
void grid_boundary<grid_type::raster, A>::set_status_at_boundary(const NodeStatus top,
                                                                 const NodeStatus right,
                                                                 const NodeStatus bottom,
                                                                 const NodeStatus left)
{
    // check symmetry of looped boundaries
    if ((top == NodeStatus::LOOPED_BOUNDARY ^ bottom == NodeStatus::LOOPED_BOUNDARY) ||
        (left == NodeStatus::LOOPED_BOUNDARY ^ right == NodeStatus::LOOPED_BOUNDARY))
    {
        throw std::invalid_argument("looped boundaries are not symmetrical");
    }

    const auto shape = this->node_status.shape();
    const auto last_row = static_cast<std::ptrdiff_t>(shape[0] - 1);
    const auto last_col = static_cast<std::ptrdiff_t>(shape[1] - 1);

    auto top_nodes = xt::view(this->node_status, 0, xt::all());
    top_nodes = utils::cast_underlying(top);

    auto right_nodes = xt::view(this->node_status, xt::all(), last_col);
    right_nodes = utils::cast_underlying(right);

    auto bottom_nodes = xt::view(this->node_status, last_row, xt::all());
    bottom_nodes = utils::cast_underlying(bottom);

    auto left_nodes = xt::view(this->node_status, xt::all(), 0);
    left_nodes = utils::cast_underlying(left);

    // maybe fix grid corners
    if (priority[right] < priority[top])
    {
        this->node_status(0, last_col) = utils::cast_underlying(top);
    }
    if (priority[bottom] < priority[right])
    {
        this->node_status(last_row, last_col) = utils::cast_underlying(right);
    }
    if (priority[left] < priority[bottom])
    {
        this->node_status(last_row, 0) = utils::cast_underlying(bottom);
    }
}


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
template<class N>
void set_node_status_grid_boundaries(N& node_status,
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


/**
 * Helper function for setting node status at each side of the grid
 * boundaries.
 *
 * @param node_status : ``[intent=inout, shape=(nrows, ncols)]``
 *     Array of node status on the grid.
 * @param border :
 *     Node status to set on all grid sides.
 */
template<class N>
void set_node_status_grid_boundaries(N& node_status,
                                     NodeStatus border)
{
    auto row_bound = xt::view(node_status, xt::keep(0, -1), xt::all());
    row_bound = border;

    auto col_bound = xt::view(node_status, xt::all(), xt::keep(0, -1));
    col_bound = border;
}


namespace detail
{


/**
 * A very weak (but efficient) way to check for looped boundary
 * conditions on rows.
 */
template<class N>
bool is_rows_looped_boundaries(N&& node_status)
{
    return node_status(0, 1) == fastscapelib::NodeStatus::LOOPED_BOUNDARY;
}


/**
 * A very weak (but efficient) way to check for looped boundary
 * conditions on cols.
 */
template<class N>
bool is_cols_looped_boundaries(N&& node_status)
{
    return node_status(1, 0) == fastscapelib::NodeStatus::LOOPED_BOUNDARY;
}


}  // namespace detail

}  // namespace fastscapelib
