/**
 * grids (channel, raster, etc.) hold and manage grid/mesh geometry,
 * topology and boundary conditions.
 */
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <vector>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xtensor.hpp"

#include "fastscapelib/utils.hpp"
#include "fastscapelib/xtensor_utils.hpp"


namespace fastscapelib
{

//*****************
//* Grid boundaries
//*****************

/**
 * Status at a grid/mesh node.
 */
enum class node_status : std::uint8_t
{
    core = 0,
    fixed_value_boundary = 1,
    fixed_gradient_boundary = 2,
    looped_boundary = 3
};

/**
 * Status at grid boundary nodes.
 */
struct boundary_status
{
    node_status left = node_status::core;    /**< Status at left edge/border node(s) */
    node_status right = node_status::core;   /**< Status at right edge/border node(s) */
    node_status top = node_status::core;     /**< Status at top border nodes */
    node_status bottom = node_status::core;  /**< Status at bottom border nodes */

    boundary_status(node_status boundaries);
    boundary_status(std::initializer_list<node_status> boundaries);
    boundary_status(const std::array<node_status, 2>& edges);
    boundary_status(const std::array<node_status, 4>& borders);
};

/**
 * @name Constructors
 */
//@{
/**
 * Set the same status for all boundary nodes.
 */
inline boundary_status::boundary_status(node_status boundaries)
    : left(boundaries), right(boundaries), top(boundaries), bottom(boundaries)
{
}

/**
 * Set status either at the left/right edges nodes on a profile grid (2-length list)
 * or at the left/right/top/bottom border nodes on a raster grid (4-length list).
 */
inline boundary_status::boundary_status(std::initializer_list<node_status> boundaries)
{
    auto bnd = boundaries.begin();
    left = *bnd++;
    right = *bnd++;
    if (boundaries.size() == 4)
    {
        top = *bnd++;
        bottom = *bnd++;
    }
}

/**
 * Set status at the left and right edge nodes of a profile grid.
 */
inline boundary_status::boundary_status(const std::array<node_status, 2>& edges)
    : left(edges[0]), right(edges[1])
{
}

/**
 * Set status at the left, right, top and bottom border nodes of a raster grid.
 */
inline boundary_status::boundary_status(const std::array<node_status, 4>& borders)
    : left(borders[0]), right(borders[1]), top(borders[2]), bottom(borders[3])
{
}
//@}

//***************
//* Grid elements
//***************

/**
 * Grid/mesh node.
 */
struct node
{
    std::size_t idx;     /**< Node index */
    node_status status;  /**< Node status */
};

/**
 * Raster grid node.
 */
struct raster_node
{
    std::size_t row;     /**< Row index */
    std::size_t col;     /**< Column index */
    node_status status;  /**< Node status */
};

/**
 * Neighbor node.
 */
struct neighbor
{
    std::size_t idx;     /**< Index of the neighbor node */
    double distance;     /**< Distance to the neighbor node */
    node_status status;  /**< Status at the neighbor node */
};

/**
 * Neighbor node (on raster grids).
 */
struct raster_neighbor
{
    std::size_t flatten_idx;  /**< Flattened index of the neighbor node */
    std::size_t row;          /**< Row index of the neighbor node */
    std::size_t col;          /**< Column index of the neighbor node */
    double distance;          /**< Distance to the neighbor node */
    node_status status;       /**< Status at the neighbor node */
};

namespace detail
{

/**
 * Add a given offset (positive, negative or zero) to a grid node index.
 *
 * This utility function is mainly to avoid signed/unsigned conversion
 * warnings. Use it carefully.
 */
inline std::size_t add_offset(std::size_t idx, std::ptrdiff_t offset)
{
    return static_cast<std::size_t>(static_cast<std::ptrdiff_t>(idx) + offset);
}

}

//****************
//* Grid interface
//****************

template <class G>
struct grid_inner_types;

/**
 * @class grid_interface
 * @brief Add a common interface to all grid types.
 *
 * This class only defines a basic interface for all grid types.
 * It does not embed any data member, this responsibility
 * is delegated to the inheriting grid classes.
 *
 * @tparam G Derived grid type.
 */
template <class G>
class grid_interface
{
public:

    using derived_grid_type = G;
    using inner_types = grid_inner_types<derived_grid_type>;
    using spacing_type = typename inner_types::spacing_type;
    //using shape_type = typename grid_type::shape_type;
    using node_status_type = typename inner_types::node_status_type;
    using neighbors_type = typename inner_types::neighbors_type;

    //shape_type shape() const noexcept;
    std::size_t size() const noexcept;
    spacing_type spacing() const noexcept;
    const node_status_type& status_at_nodes() const;

    const neighbors_type& neighbors(std::size_t idx) const;

protected:
    grid_interface() = default;
    ~grid_interface() = default;

    const derived_grid_type& derived_grid() const & noexcept;
};

template <class G>
inline auto grid_interface<G>::derived_grid() const & noexcept -> const derived_grid_type&
{
    return *static_cast<const derived_grid_type*>(this);
}

/**
 * @name Grid properties
 */
//@{
/**
 * Returns the total number of grid nodes.
 */
template <class G>
inline std::size_t grid_interface<G>::size() const noexcept
{
    return derived_grid().m_size;
}

/**
 * Returns the (uniform) spacing between two adjacent grid nodes.
 */
template <class G>
inline auto grid_interface<G>::spacing() const noexcept -> spacing_type
{
    return derived_grid().m_spacing;
}

/**
 * Returns a const reference to the array of status at grid nodes.
 */
template <class G>
inline auto grid_interface<G>::status_at_nodes() const -> const node_status_type&
{
    return derived_grid().m_status_at_nodes;
}
//@}

/**
 * @name Iterators
 */
/**
 * Iterate over the neighbors of a given grid node.
 *
 * Follows looped boundary conditions, if any.
 *
 * @param idx Index of the grid node.
 * @return Reference to the vector of the neighbors of that grid node.
 */
template <class G>
inline auto grid_interface<G>::neighbors(std::size_t idx) const -> const neighbors_type&
{
    return derived_grid().neighbors_impl(idx);
}
//@}


//*******************
//* Profile grid (1D)
//*******************

template <class XT>
class profile_grid_xt;

template <class XT>
struct grid_inner_types<profile_grid_xt<XT>>
{
    static constexpr bool is_structured = true;
    static constexpr bool is_uniform = true;

    using spacing_type = double;
    using node_status_type = xt_container_t<XT, node_status, 1>;
    using neighbors_type = std::vector<neighbor>;
};

/**
 * @class profile_grid_xt
 * @brief 1-dimensional uniform grid.
 *
 * Used for modeling single channel or hillslope profiles.
 *
 * @tparam XT xtensor container selector for data array members.
 */
template <class XT>
class profile_grid_xt : public grid_interface<profile_grid_xt<XT>>
{
public:

    using self_type = profile_grid_xt<XT>;
    using base_type = grid_interface<self_type>;

    using spacing_type = typename base_type::spacing_type;
    using node_status_type = typename base_type::node_status_type;
    using neighbors_type = typename base_type::neighbors_type;

    profile_grid_xt(std::size_t size,
                    double spacing,
                    const boundary_status& status_at_bounds,
                    const std::vector<node>& status_at_nodes = {});

private:

    std::size_t m_size;
    double m_spacing;

    node_status_type m_status_at_nodes;
    boundary_status m_status_at_bounds;
    bool has_looped_edges = false;
    void set_status_at_nodes(const std::vector<node>& status_at_nodes);

    static const std::array<std::ptrdiff_t, 3> offsets;
    std::vector<neighbors_type> m_all_neighbors;
    void precompute_neighbors();

    const neighbors_type& neighbors_impl(std::size_t idx) const;

    friend class grid_interface<self_type>;
};

template <class XT>
constexpr std::array<std::ptrdiff_t, 3> profile_grid_xt<XT>::offsets = { {0, -1, 1} };

/**
 * @name Constructors
 */
//@{
/**
 * Creates a new profile grid.
 *
 * @param size Total number of grid nodes.
 * @param spacing Distance between two adjacent grid nodes.
 * @param status_at_bounds Status at boundary nodes (left & right grid edges).
 * @param status_at_nodes Manually define the status at any node on the grid.
 */
template <class XT>
profile_grid_xt<XT>::profile_grid_xt(std::size_t size,
                                     double spacing,
                                     const boundary_status& status_at_bounds,
                                     const std::vector<node>& status_at_nodes)
    : base_type(), m_size(size), m_spacing(spacing), m_status_at_bounds(status_at_bounds)
{
    set_status_at_nodes(status_at_nodes);
    precompute_neighbors();
}
//@}

template <class XT>
void profile_grid_xt<XT>::set_status_at_nodes(const std::vector<node>& status_at_nodes)
{
    std::array<std::size_t, 1> shape {m_size};
    m_status_at_nodes.resize(shape);
    m_status_at_nodes.fill(node_status::core);

    m_status_at_nodes[0] = m_status_at_bounds.left;
    m_status_at_nodes[m_size-1] = m_status_at_bounds.right;

    for (const node& inode : status_at_nodes)
    {
        m_status_at_nodes(inode.idx) = inode.status;
    }

    bool left_looped = m_status_at_nodes[0] == node_status::looped_boundary;
    bool right_looped = m_status_at_nodes[m_size-1] == node_status::looped_boundary;

    if (left_looped ^ right_looped)
    {
        throw std::invalid_argument("inconsistent looped boundary status at grid edges");
    }
    else if (left_looped && right_looped)
    {
        has_looped_edges = true;
    }
}

template <class XT>
void profile_grid_xt<XT>::precompute_neighbors()
{
    m_all_neighbors.resize(m_size);

    for (std::size_t idx=1; idx<m_size-1; ++idx)
    {
        for (std::size_t k=1; k<3; ++k)
        {
            std::size_t nb_idx = detail::add_offset(idx, offsets[k]);
            neighbor nb {nb_idx, m_spacing, m_status_at_nodes[nb_idx]};
            m_all_neighbors[idx].push_back(nb);
        }
    }

    m_all_neighbors[0].push_back(
        {1, m_spacing, m_status_at_nodes[1]});
    m_all_neighbors[m_size-1].push_back(
        {m_size-2, m_spacing, m_status_at_nodes[m_size-2]});

    if (has_looped_edges)
    {
        m_all_neighbors[0].insert(
            m_all_neighbors[0].begin(),
            {m_size-1, m_spacing, m_status_at_nodes[m_size-1]});
        m_all_neighbors[m_size-1].push_back(
            {0, m_spacing, m_status_at_nodes[0]});
    }
}

template <class XT>
inline auto profile_grid_xt<XT>::neighbors_impl(std::size_t idx) const
    -> const neighbors_type&
{
    return m_all_neighbors[idx];
}

/**
 * @typedef profile_grid
 * Alias template on profile_grid_xt with ``xt::xtensor`` used
 * as array container type for data members.
 *
 * This is mainly for convenience when using in C++ applications.
 *
 */
using profile_grid = profile_grid_xt<xtensor_selector>;


//******************
//* Raster grid (2D)
//******************

enum class raster_connect : std::size_t
{
   king = 4,
   queen = 8
};

} // namespace fastscapelib
