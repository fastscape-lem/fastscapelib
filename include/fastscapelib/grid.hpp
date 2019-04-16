/**
 * grids (channel, raster, etc.) hold and manage grid/mesh geometry,
 * topology and boundary conditions.
 */
#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <type_traits>
#include <map>
#include <vector>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xindex_view.hpp"

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


namespace detail
{

inline bool node_status_cmp(node_status a, node_status b)
{
    static std::map<node_status, int> priority {
        {node_status::core, 0},
        {node_status::looped_boundary, 1},
        {node_status::fixed_gradient_boundary, 2},
        {node_status::fixed_value_boundary, 3}
    };

    return priority[a] < priority[b];
}

}  // namespace detail


/**
 * Status at grid boundary nodes.
 */
class boundary_status
{
public:

    node_status left = node_status::core;    /**< Status at left edge/border node(s) */
    node_status right = node_status::core;   /**< Status at right edge/border node(s) */
    node_status top = node_status::core;     /**< Status at top border nodes */
    node_status bottom = node_status::core;  /**< Status at bottom border nodes */

    boundary_status(node_status boundaries);
    boundary_status(std::initializer_list<node_status> boundaries);
    boundary_status(const std::array<node_status, 2>& edges);
    boundary_status(const std::array<node_status, 4>& borders);

    bool is_horizontal_looped();
    bool is_vertical_looped();

private:

    bool is_looped(node_status status);
    void check_looped_symmetrical();
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
    check_looped_symmetrical();
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

    check_looped_symmetrical();
}

/**
 * Set status at the left and right edge nodes of a profile grid.
 */
inline boundary_status::boundary_status(const std::array<node_status, 2>& edges)
    : left(edges[0]), right(edges[1])
{
    check_looped_symmetrical();
}

/**
 * Set status at the left, right, top and bottom border nodes of a raster grid.
 */
inline boundary_status::boundary_status(const std::array<node_status, 4>& borders)
    : left(borders[0]), right(borders[1]), top(borders[2]), bottom(borders[3])
{
    check_looped_symmetrical();
}
//@}

inline bool boundary_status::is_looped(node_status status)
{
    return status == node_status::looped_boundary;
}

/**
 * Return true if periodic (looped) conditions are set for ``left`` and ``right``.
 */
inline bool boundary_status::is_horizontal_looped()
{
    return is_looped(left) && is_looped(right);
}

/**
 * Return true if periodic (looped) conditions are set for ``top`` and ``bottom``.
 */
inline bool boundary_status::is_vertical_looped()
{
    return is_looped(top) && is_looped(bottom);
}

inline void boundary_status::check_looped_symmetrical()
{
    if (is_looped(left) ^ is_looped(right) || is_looped(top) ^ is_looped(bottom))
    {
        throw std::invalid_argument("looped boundaries are not symmetrical");
    }
}

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

    using size_type = typename inner_types::size_type;
    using shape_type = typename inner_types::shape_type;
    using spacing_type = typename inner_types::spacing_type;
    using neighbors_type = typename inner_types::neighbors_type;

    using spacing_t = std::conditional_t<std::is_arithmetic<spacing_type>::value,
                                         spacing_type,
                                         const spacing_type&>;

    using node_status_type = typename inner_types::node_status_type;

    //shape_type shape() const noexcept;
    size_type size() const noexcept;
    spacing_t spacing() const noexcept;
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
inline auto grid_interface<G>::size() const noexcept -> size_type
{
    return derived_grid().m_size;
}

/**
 * Returns the (uniform) spacing between two adjacent grid nodes.
 *
 * Returns a copy of the value for 1-d grids or a constant reference otherwise.
 */
template <class G>

inline auto grid_interface<G>::spacing() const noexcept -> spacing_t
{
    return derived_grid().m_spacing;
}

/**
 * Returns a constant reference to the array of status at grid nodes.
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
    using xt_selector = XT;

    struct xt_ndims
    {
        static constexpr std::size_t value = 1;
    };

    using xt_type = xt_container_t<xt_selector, int, xt_ndims::value>;

    using size_type = typename xt_type::size_type;
    using shape_type = typename xt_type::shape_type;

    struct is_structured
    {
        static constexpr bool value = true;
    };

    struct is_uniform
    {
        static constexpr bool value = true;
    };

    using spacing_type = double;
    using node_status_type = xt_container_t<xt_selector, node_status, xt_ndims::value>;
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
    using inner_types = grid_inner_types<self_type>;

    using xt_selector = typename inner_types::xt_selector;
    using xt_ndims = typename inner_types::xt_ndims;

    using size_type = typename inner_types::size_type;
    using shape_type = typename inner_types::shape_type;

    using spacing_type = typename inner_types::spacing_type;
    using neighbors_type = typename inner_types::neighbors_type;

    using node_status_type = typename base_type::node_status_type;

    profile_grid_xt(size_type size,
                    spacing_type spacing,
                    const boundary_status& status_at_bounds,
                    const std::vector<node>& status_at_nodes = {});

private:

    shape_type m_shape;
    size_type m_size;
    spacing_type m_spacing;

    node_status_type m_status_at_nodes;
    boundary_status m_status_at_bounds;
    void set_status_at_nodes(const std::vector<node>& status_at_nodes);

    static constexpr std::array<std::ptrdiff_t, 3> offsets { {0, -1, 1} };
    std::vector<neighbors_type> m_all_neighbors;
    void precompute_neighbors();

    const neighbors_type& neighbors_impl(std::size_t idx) const;

    friend class grid_interface<self_type>;
};

template <class XT>
constexpr std::array<std::ptrdiff_t, 3> profile_grid_xt<XT>::offsets;

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
profile_grid_xt<XT>::profile_grid_xt(size_type size,
                                     spacing_type spacing,
                                     const boundary_status& status_at_bounds,
                                     const std::vector<node>& status_at_nodes)
    : base_type(), m_size(size), m_spacing(spacing), m_status_at_bounds(status_at_bounds)
{
    m_shape = {static_cast<typename shape_type::value_type>(m_size)};
    set_status_at_nodes(status_at_nodes);
    precompute_neighbors();
}
//@}

template <class XT>
void profile_grid_xt<XT>::set_status_at_nodes(const std::vector<node>& status_at_nodes)
{
    node_status_type temp_status_at_nodes(m_shape, node_status::core);

    temp_status_at_nodes(0) = m_status_at_bounds.left;
    temp_status_at_nodes(m_size-1) = m_status_at_bounds.right;

    for (const node& n : status_at_nodes)
    {
        if (n.status == node_status::looped_boundary)
        {
            throw std::invalid_argument("node_status::looped_boundary is not allowed in "
                                        "'status_at_nodes' "
                                        "(use 'status_at_bounds' instead)");
        }
        else if (temp_status_at_nodes(n.idx) == node_status::looped_boundary)
        {
            throw std::invalid_argument("cannot overwrite the status of a "
                                        "looped boundary node");
        }

        temp_status_at_nodes(n.idx) = n.status;
    }

    m_status_at_nodes = temp_status_at_nodes;
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

    if (m_status_at_bounds.is_horizontal_looped())
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

/**
 * Raster grid node connectivity.
 */
enum class raster_connect
{
    rook = 0,   /**< 4-connectivity (vertical or horizontal) */
    queen       /**< 8-connectivity (including diagonals) */
    // bishop   /**< 4-connectivity (only diagonals) */
};

template <class XT>
class raster_grid_xt;

template <class XT>
struct grid_inner_types<raster_grid_xt<XT>>
{
    using xt_selector = XT;

    struct xt_ndims
    {
        static constexpr std::size_t value = 2;
    };

    using xt_type = xt_container_t<xt_selector, int, xt_ndims::value>;

    using size_type = typename xt_type::size_type;
    using shape_type = typename xt_type::shape_type;

    struct is_structured
    {
        static constexpr bool value = true;
    };

    struct is_uniform
    {
        static constexpr bool value = true;
    };

    using spacing_type = std::array<double, 2>;
    using node_status_type = xt_container_t<xt_selector, node_status, xt_ndims::value>;
    using neighbors_type = std::vector<raster_neighbor>;
};


/**
 * @class raster_grid_xt
 * @brief 2-dimensional uniform (raster) grid.
 *
 * @tparam XT xtensor container selector for data array members.
 */
template <class XT>
class raster_grid_xt : public grid_interface<raster_grid_xt<XT>>
{
public:

    using self_type = profile_grid_xt<XT>;
    using base_type = grid_interface<self_type>;
    using inner_types = grid_inner_types<self_type>;

    using xt_selector = typename inner_types::xt_selector;
    using xt_ndims = typename inner_types::xt_ndims;

    using size_type = typename inner_types::size_type;
    using shape_type = typename inner_types::shape_type;

    using spacing_type = typename inner_types::spacing_type;
    using neighbors_type = typename inner_types::neighbors_type;
    using neighbor_indices_type = xt::xtensor<std::array<size_type, 2>, 1>;

    using node_status_type = typename inner_types::node_status_type;

    raster_grid_xt(const shape_type& shape,
                   const spacing_type& spacing,
                   const boundary_status& status_at_bounds,
                   const std::vector<raster_node>& status_at_nodes = {});

    template <raster_connect RC = raster_connect::queen>
    inline const neighbor_indices_type neighbor_indices(size_type row,
                                                        size_type col) const noexcept;

    template <class E, raster_connect RC = raster_connect::queen>
    inline auto neighbor_view(E&& field, size_type row, size_type col) noexcept;

private:

    shape_type m_shape;
    size_type m_size;
    spacing_type m_spacing;

    node_status_type m_status_at_nodes;
    boundary_status m_status_at_bounds;
    void set_status_at_nodes(const std::vector<raster_node>& status_at_nodes);

    struct corner_node
    {
        size_type row;
        size_type col;
        node_status row_border;
        node_status col_border;
    };

    std::array<std::vector<std::uint8_t>, 2> m_gcode_components;
    void build_gcode();
    std::uint8_t gcode(size_type row, size_type col) const noexcept;

    using offset_list = xt::xtensor<std::array<std::ptrdiff_t, 2>, 1>;

    template <raster_connect RC, std::enable_if_t<RC == raster_connect::queen, int> = 0>
    offset_list make_offsets_list(std::ptrdiff_t up, std::ptrdiff_t down,
                                   std::ptrdiff_t left, std::ptrdiff_t right) const;

    template <raster_connect RC, std::enable_if_t<RC == raster_connect::rook, int> = 0>
    offset_list make_offsets_list(std::ptrdiff_t up, std::ptrdiff_t down,
                                   std::ptrdiff_t left, std::ptrdiff_t right) const;

    using neighbor_offsets_type = std::array<offset_list, 9>;

    neighbor_offsets_type m_neighbor_offsets_queen;
    neighbor_offsets_type m_neighbor_offsets_rook;

    template <raster_connect RC>
    neighbor_offsets_type build_neighbor_offsets();

    template <raster_connect RC, std::enable_if_t<RC == raster_connect::queen, int> = 0>
    inline const offset_list& neighbor_offsets(size_type row, size_type col) const noexcept;

    template <raster_connect RC, std::enable_if_t<RC == raster_connect::rook, int> = 0>
    inline const offset_list& neighbor_offsets(size_type row, size_type col) const noexcept;

    const neighbors_type& neighbors_impl(std::size_t row, std::size_t col) const;
    const neighbors_type& neighbors_impl(std::size_t idx) const;
};

/**
 * @name Constructors
 */
//@{
/**
 * Creates a new raster grid.
 *
 * @param shape Shape of the grid (number of rows and cols).
 * @param spacing Distance between two adjacent rows / cols.
 * @param status_at_bounds Status at boundary nodes (grid borders).
 * @param status_at_nodes Manually define the status at any node on the grid.
 */
template <class XT>
raster_grid_xt<XT>::raster_grid_xt(const shape_type& shape,
                                   const spacing_type& spacing,
                                   const boundary_status& status_at_bounds,
                                   const std::vector<raster_node>& status_at_nodes)
    : m_shape(shape), m_spacing(spacing), m_status_at_bounds(status_at_bounds)
{
    m_size = shape[0] * shape[1];

    build_gcode();

    set_status_at_nodes(status_at_nodes);

    m_neighbor_offsets_queen = build_neighbor_offsets<raster_connect::queen>();
    m_neighbor_offsets_rook = build_neighbor_offsets<raster_connect::rook>();
}
//@}

template <class XT>
void raster_grid_xt<XT>::set_status_at_nodes(const std::vector<raster_node>& status_at_nodes)
{
    node_status_type temp_status_at_nodes(m_shape, node_status::core);

    // set border nodes
    auto left = xt::view(temp_status_at_nodes, xt::all(), 0);
    left = m_status_at_bounds.left;

    auto right = xt::view(temp_status_at_nodes, xt::all(), xt::keep(-1));
    right = m_status_at_bounds.right;

    auto top = xt::view(temp_status_at_nodes, 0, xt::all());
    top = m_status_at_bounds.top;

    auto bottom = xt::view(temp_status_at_nodes, xt::keep(-1), xt::all());
    bottom = m_status_at_bounds.bottom;

    // set corner nodes
    std::vector<corner_node> corners = {
        {0, 0, m_status_at_bounds.top, m_status_at_bounds.left},
        {0, m_shape[1]-1, m_status_at_bounds.top, m_status_at_bounds.right},
        {m_shape[0]-1, 0, m_status_at_bounds.bottom, m_status_at_bounds.left},
        {m_shape[0]-1, m_shape[1]-1, m_status_at_bounds.bottom, m_status_at_bounds.right}
    };

    for (const auto& c : corners)
    {
        node_status cs = std::max(c.row_border, c.col_border, detail::node_status_cmp);
        temp_status_at_nodes(c.row, c.col) = cs;
    }

    // set user-defined nodes
    for (const raster_node& n : status_at_nodes)
    {
        if (n.status == node_status::looped_boundary)
        {
            throw std::invalid_argument("node_status::looped_boundary is not allowed in "
                                        "'status_at_nodes' "
                                        "(use 'status_at_bounds' instead)");
        }
        else if (temp_status_at_nodes(n.row, n.col) == node_status::looped_boundary)
        {
            throw std::invalid_argument("cannot overwrite the status of a "
                                        "looped boundary node");
        }

        temp_status_at_nodes(n.row, n.col) = n.status;
    }

    m_status_at_nodes = temp_status_at_nodes;
}

/**
 * Pre-store for each (row, col) dimension 1-d arrays
 * that will be used to get the characteristic location of a node
 * on the grid.
 */
template <class XT>
void raster_grid_xt<XT>::build_gcode()
{
    for (std::uint8_t dim=0; dim<2; ++dim)
    {
        std::uint8_t fill_value = dim * 2 + 1;
        std::vector<std::uint8_t> gcode_component(m_shape[dim], fill_value);

        gcode_component[0] = 0;
        gcode_component[m_shape[dim]-1] = fill_value * 2;

        m_gcode_components[dim] = gcode_component;
    }
}

/**
 * Given row and col indexes, return a code in the range [0,8], which
 * corresponds to one of the following characteristic locations on the
 * grid (i.e., inner/border/corner):
 *
 *   0 -- 3 -- 6
 *   |         |
 *   1    4    7
 *   |         |
 *   2 -- 5 -- 8
 */
template <class XT>
inline std::uint8_t raster_grid_xt<XT>::gcode(size_type row, size_type col) const noexcept
{
    return m_gcode_components[0][row] + m_gcode_components[1][col];
}

template <class XT>
template <raster_connect RC, std::enable_if_t<RC == raster_connect::queen, int>>
auto raster_grid_xt<XT>::make_offsets_list(std::ptrdiff_t up,
                                           std::ptrdiff_t down,
                                           std::ptrdiff_t left,
                                           std::ptrdiff_t right)
    const -> offset_list
{
    xt::xtensor<bool, 1> mask {
        true,
        up != 0 || left != 0, up != 0, up != 0 || right != 0,
        left != 0, right != 0,
        down != 0 || left != 0, down != 0, down != 0 || right != 0
    };

    offset_list full_offsets {
        {0, 0},
        {up, left}, {up, 0}, {up, right},
        {0, left}, {0, right},
        {down, left}, {down, 0}, {down, right}
    };

    return xt::filter(full_offsets, mask);
}

template <class XT>
template <raster_connect RC, std::enable_if_t<RC == raster_connect::rook, int>>
auto raster_grid_xt<XT>::make_offsets_list(std::ptrdiff_t up,
                                           std::ptrdiff_t down,
                                           std::ptrdiff_t left,
                                           std::ptrdiff_t right)
    const -> offset_list
{
    xt::xtensor<bool, 1> mask {
        true, up != 0, left != 0, right != 0, down != 0
    };

    offset_list full_offsets {
        {0, 0}, {up, 0}, {0, left}, {0, right}, {down, 0}
    };

    return xt::filter(full_offsets, mask);
}

/**
 * Pre-store the (row, col) index offsets of grid node neighbors for
 * each of the 9-characteristic locations on the grid.
 *
 * Those offsets follow any looped conditions given on grid
 * boundaries. First offset is zero (current node) and further offsets
 * are stored following a row-major layout.
 */
template <class XT>
template <raster_connect RC>
auto raster_grid_xt<XT>::build_neighbor_offsets() -> neighbor_offsets_type
{
    auto dr = m_shape[0] - 1;
    auto dc = m_shape[1] - 1;

    if (!m_status_at_bounds.is_vertical_looped())
    {
        dr = 0;
    }
    if (!m_status_at_bounds.is_horizontal_looped())
    {
        dc = 0;
    }

    return {
        make_offsets_list<RC>(dr, 1, dc, 1),      // top-left corner
        make_offsets_list<RC>(-1, 1, dc, 1),      // left border
        make_offsets_list<RC>(-1, -dr, dc, 1),    // bottom-left corner
        make_offsets_list<RC>(dr, 1, -1, 1),      // top border
        make_offsets_list<RC>(-1, 1, -1, 1),      // inner
        make_offsets_list<RC>(-1, -dr, -1, 1),    // bottom border
        make_offsets_list<RC>(dr, 1, -1, -dc),    // top-right corner
        make_offsets_list<RC>(-1, 1, -1, -dc),    // right border
        make_offsets_list<RC>(-1, -dr, -1, -dc),  // bottom-right corner
    };
}

template <class XT>
template <raster_connect RC, std::enable_if_t<RC == raster_connect::queen, int>>
inline auto raster_grid_xt<XT>::neighbor_offsets(size_type row,
                                                 size_type col) const noexcept
    -> const offset_list&
{
    return m_neighbor_offsets_queen(gcode(row, col));
}

template <class XT>
template <raster_connect RC, std::enable_if_t<RC == raster_connect::rook, int>>
inline auto raster_grid_xt<XT>::neighbor_offsets(size_type row,
                                                 size_type col) const noexcept
    -> const offset_list&
{
    return m_neighbor_offsets_rook(gcode(row, col));
}

template <class XT>
template <raster_connect RC>
inline auto raster_grid_xt<XT>::neighbor_indices(size_type row,
                                                 size_type col) const noexcept
    -> const neighbor_indices_type
{
    auto offsets = neighbor_offsets<RC>(row, col);

}

template <class XT>
template <class E, raster_connect RC>
inline auto raster_grid_xt<XT>::neighbor_view(E&& field,
                                              size_type row,
                                              size_type col) noexcept
{
    return xt::index_view(std::forward<E>(field), neighbor_indices<RC>(row, col));
}

/**
 * @typedef raster_grid
 * Alias template on raster_grid_xt with ``xt::xtensor`` used
 * as array container type for data members.
 *
 * This is mainly for convenience when using in C++ applications.
 *
 */
using raster_grid = raster_grid_xt<xtensor_selector>;


} // namespace fastscapelib
