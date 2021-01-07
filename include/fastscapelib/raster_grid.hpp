/**
 * A raster grid represents a 2D uniform grid.
 */
#pragma once

#include "fastscapelib/structured_grid.hpp"
#include "fastscapelib/profile_grid.hpp"


namespace fastscapelib
{
    //***************
    //* Grid elements
    //***************

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
     * Neighbor node (on raster grids).
     */
    struct raster_neighbor
    {
        std::size_t flatten_idx;  /**< Flattened index of the neighbor node */
        std::size_t row;          /**< Row index of the neighbor node */
        std::size_t col;          /**< Column index of the neighbor node */
        double distance;          /**< Distance to the neighbor node */
        node_status status;       /**< Status at the neighbor node */

        bool operator==(const raster_neighbor& rhs) const
        {
            return (flatten_idx == rhs.flatten_idx) && (row == rhs.row) && (col == rhs.col) && 
                    (distance == rhs.distance) && (status == rhs.status);
        }
    };

    /**
     * Boundary status (on raster grids).
     */
    class raster_boundary_status: public profile_boundary_status
    {
    public:

        node_status top = node_status::core;     /**< Status at top border nodes */
        node_status bottom = node_status::core;  /**< Status at bottom border nodes */

        raster_boundary_status(node_status status);
        raster_boundary_status(const std::array<node_status, 4>& status);

        bool is_vertical_looped();

    protected:

        void check_looped_symmetrical();        
    };

    /**
     * @name Constructors
     */
    //@{
    /**
     * Set the same status for all boundary nodes.
     */
    inline raster_boundary_status::raster_boundary_status(node_status status)
        : profile_boundary_status(status), top(status), bottom(status)
    {
        check_looped_symmetrical();
    }

    /**
     * Set status at the left, right, top and bottom border nodes of a raster grid.
     */
    inline raster_boundary_status::raster_boundary_status(const std::array<node_status, 4>& status)
        : profile_boundary_status(status[0], status[1]), top(status[2]), bottom(status[3])
    {
        check_looped_symmetrical();
    }
    //@}

    /**
     * Return true if periodic (looped) conditions are set for ``top`` and ``bottom``.
     */
    inline bool raster_boundary_status::is_vertical_looped()
    {
        return is_looped(top) && is_looped(bottom);
    }

    inline void raster_boundary_status::check_looped_symmetrical()
    {
        if (is_looped(left) ^ is_looped(right) || is_looped(top) ^ is_looped(bottom))
        {
            throw std::invalid_argument("looped boundaries are not symmetrical");
        }
    }

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

        using boundary_status_type = raster_boundary_status;
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
    class raster_grid_xt : public structured_grid<raster_grid_xt<XT>>
    {
    public:

        using self_type = raster_grid_xt<XT>;
        using base_type = structured_grid<self_type>;
        using inner_types = grid_inner_types<self_type>;

        using xt_selector = typename inner_types::xt_selector;
        using xt_ndims = typename inner_types::xt_ndims;

        using size_type = typename inner_types::size_type;
        using shape_type = typename inner_types::shape_type;

        using spacing_type = typename inner_types::spacing_type;
        using neighbors_type = typename inner_types::neighbors_type;
        using neighbor_indices_type = xt::xtensor<
            xt::xtensor_fixed<size_type, xt::xshape<2>>, 1>;
        using distance_list = xt::xtensor<double, 1>;

        using boundary_status_type = typename inner_types::boundary_status_type;
        using node_status_type = typename inner_types::node_status_type;

        template <class S>
        raster_grid_xt(S&& shape,
                    const spacing_type& spacing,
                    const boundary_status_type& status_at_bounds,
                    const std::vector<raster_node>& status_at_nodes = {});

        template <raster_connect RC = raster_connect::queen>
        inline neighbor_indices_type neighbor_indices(size_type row,
                                                            size_type col) const noexcept;

        template <raster_connect RC = raster_connect::queen>
        inline const distance_list& neighbor_distances(size_type row,
                                                    size_type col) const noexcept;

        template <raster_connect RC = raster_connect::queen, class E>
        inline auto neighbor_view(E&& field, size_type row, size_type col) const noexcept;

        std::uint8_t gcode(size_type row, size_type col) const noexcept;

        using offset_list = xt::xtensor<xt::xtensor_fixed<std::ptrdiff_t, xt::xshape<2>>, 1>;

        template <raster_connect RC, std::enable_if_t<RC == raster_connect::queen, int> = 0>
        inline const offset_list& neighbor_offsets(std::uint8_t code) const noexcept;

        template <raster_connect RC, std::enable_if_t<RC == raster_connect::rook, int> = 0>
        inline const offset_list& neighbor_offsets(std::uint8_t code) const noexcept;

    private:

        shape_type m_shape;
        size_type m_size;
        spacing_type m_spacing;

        node_status_type m_status_at_nodes;
        boundary_status_type m_status_at_bounds;
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


        using neighbor_distances_type = std::array<distance_list, 9>;

        neighbor_distances_type m_neighbor_distances_queen;
        neighbor_distances_type m_neighbor_distances_rook;

        template <raster_connect RC>
        neighbor_distances_type build_neighbor_distances();

        template <raster_connect RC, std::enable_if_t<RC == raster_connect::queen, int> = 0>
        inline const distance_list& neighbor_distances_impl(size_type row,
                                                            size_type col) const noexcept;

        template <raster_connect RC, std::enable_if_t<RC == raster_connect::rook, int> = 0>
        inline const distance_list& neighbor_distances_impl(size_type row,
                                                            size_type col) const noexcept;

        const neighbors_type& neighbors_impl(std::size_t row, std::size_t col) const;
        const neighbors_type& neighbors_impl(std::size_t idx) const;

        friend class structured_grid<self_type>;
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
    template <class S>
    raster_grid_xt<XT>::raster_grid_xt(S&& shape,
                                    const spacing_type& spacing,
                                    const boundary_status_type& status_at_bounds,
                                    const std::vector<raster_node>& status_at_nodes)
        : m_shape(shape), m_spacing(spacing), m_status_at_bounds(status_at_bounds)
    {
        m_size = shape[0] * shape[1];

        build_gcode();

        set_status_at_nodes(status_at_nodes);

        m_neighbor_offsets_queen = build_neighbor_offsets<raster_connect::queen>();
        m_neighbor_offsets_rook = build_neighbor_offsets<raster_connect::rook>();

        m_neighbor_distances_queen = build_neighbor_distances<raster_connect::queen>();
        m_neighbor_distances_rook = build_neighbor_distances<raster_connect::rook>();
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
            else if (temp_status_at_nodes.at(n.row, n.col) == node_status::looped_boundary)
            {
                throw std::invalid_argument("cannot overwrite the status of a "
                                            "looped boundary node");
            }

            temp_status_at_nodes.at(n.row, n.col) = n.status;
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
            auto fill_value = static_cast<std::uint8_t>(3 - dim * 2);
            std::vector<std::uint8_t> gcode_component(m_shape[dim], fill_value);

            gcode_component[0] = 0;
            gcode_component[m_shape[dim]-1] = static_cast<std::uint8_t>(fill_value * 2);

            m_gcode_components[dim] = gcode_component;
        }
    }

    /**
     * Given row and col indexes, return a code in the range [0,8], which
     * corresponds to one of the following characteristic locations on the
     * grid (i.e., inner/border/corner). Use a row-major layout:
     *
     *   0 -- 1 -- 2
     *   |         |
     *   3    4    5
     *   |         |
     *   6 -- 7 -- 8
     */
    template <class XT>
    inline std::uint8_t raster_grid_xt<XT>::gcode(size_type row, size_type col) const noexcept
    {
        return static_cast<std::uint8_t>(m_gcode_components[0][row] + m_gcode_components[1][col]);
    }


    /**
     * Given the possible moves (up/down/left/right) depending on grid periodicity,
     * return a list of offsets corresponding to "Queen" neighbors relative locations.
     * Offsets are returned in row-major layout when accessible (else dropped).
     *
     * Examples:
     * 
     *         All          w/o left       w/o top  
     *      0   1   2        0   1
     *        \ | /          | / 
     *      3 - N - 4        N - 2        0 - N - 1
     *        / | \          | \            / | \  
     *      5   6   7        3   4        2   3   4
     */
    template <class XT>
    template <raster_connect RC, std::enable_if_t<RC == raster_connect::queen, int>>
    auto raster_grid_xt<XT>::make_offsets_list(std::ptrdiff_t up,
                                               std::ptrdiff_t down,
                                               std::ptrdiff_t left,
                                               std::ptrdiff_t right)
        const -> offset_list
    {
        xt::xtensor<bool, 1> mask {
            (up != 0 && left != 0), up != 0, (up != 0 && right != 0),
            left != 0, right != 0,
            (down != 0 && left != 0), down != 0, (down != 0 && right != 0)
        };

        offset_list full_offsets {
            {up, left}, {up, 0}, {up, right},
            {0, left}, {0, right},
            {down, left}, {down, 0}, {down, right}
        };

        return xt::filter(full_offsets, mask);
    }

    /**
     * Given the possible moves (up/down/left/right) depending on grid periodicity,
     * return a list of offsets corresponding to "Rook" neighbors relative locations.
     * Offsets are returned in row-major layout when accessible (else dropped).
     *
     * Examples:
     * 
     *       All         w/o left        w/o top         etc.
     *        0            0 
     *        |            |     
     *   1 -- N -- 2       N -- 1      0 -- N -- 1
     *        |            |                |
     *        3            2                2
     */
    template <class XT>
    template <raster_connect RC, std::enable_if_t<RC == raster_connect::rook, int>>
    auto raster_grid_xt<XT>::make_offsets_list(std::ptrdiff_t up,
                                               std::ptrdiff_t down,
                                               std::ptrdiff_t left,
                                               std::ptrdiff_t right)
        const -> offset_list
    {
        xt::xtensor<bool, 1> mask {
            up != 0, left != 0, right != 0, down != 0
        };

        offset_list full_offsets {
            {up, 0}, {0, left}, {0, right}, {down, 0}
        };

        return xt::filter(full_offsets, mask);
    }

    /**
     * Pre-store the (row, col) index offsets of grid node neighbors for
     * each of the 9-characteristic locations on the grid.
     *
     * Those offsets take into account looped boundary conditions (if any).
     * The order of the returned offsets corresponds to the row-major layout.
     */
    template <class XT>
    template <raster_connect RC>
    auto raster_grid_xt<XT>::build_neighbor_offsets() -> neighbor_offsets_type
    {
        auto dr = static_cast<std::ptrdiff_t>(m_shape[0] - 1);
        auto dc = static_cast<std::ptrdiff_t>(m_shape[1] - 1);

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
            make_offsets_list<RC>(dr, 1, -1, 1),      // top border
            make_offsets_list<RC>(dr, 1, -1, -dc),    // top-right corner
            make_offsets_list<RC>(-1, 1, dc, 1),      // left border
            make_offsets_list<RC>(-1, 1, -1, 1),      // inner
            make_offsets_list<RC>(-1, 1, -1, -dc),    // right border
            make_offsets_list<RC>(-1, -dr, dc, 1),    // bottom-left corner
            make_offsets_list<RC>(-1, -dr, -1, 1),    // bottom border
            make_offsets_list<RC>(-1, -dr, -1, -dc),  // bottom-right corner
        };
    }

    template <class XT>
    template <raster_connect RC, std::enable_if_t<RC == raster_connect::queen, int>>
    inline auto raster_grid_xt<XT>::neighbor_offsets(std::uint8_t code) const noexcept
        -> const offset_list&
    {
        return m_neighbor_offsets_queen[code];
    }

    template <class XT>
    template <raster_connect RC, std::enable_if_t<RC == raster_connect::rook, int>>
    inline auto raster_grid_xt<XT>::neighbor_offsets(std::uint8_t code) const noexcept
        -> const offset_list&
    {
        return m_neighbor_offsets_rook[code];
    }

    template <class XT>
    template <raster_connect RC>
    auto raster_grid_xt<XT>::build_neighbor_distances() -> neighbor_distances_type
    {
        neighbor_distances_type nb_distances;
        xt::xtensor_fixed<double, xt::xshape<2>> xspacing {m_spacing[0], m_spacing[1]};

        auto to_dist = [&xspacing](auto&& offset) -> double {
                        auto drc = xt::where(xt::equal(offset, 0), 0., 1.) * xspacing;
                        return std::sqrt(xt::sum(xt::square(drc))(0));
                    };

        for(std::uint8_t k=0; k<9; ++k)
        {
            auto offsets = neighbor_offsets<RC>(k);
            auto distances = distance_list(offsets.shape());

            std::transform(offsets.cbegin(), offsets.cend(), distances.begin(), to_dist);

            nb_distances[k] = distances;
        }

        return nb_distances;
    }

    template <class XT>
    template <raster_connect RC, std::enable_if_t<RC == raster_connect::queen, int>>
    inline auto raster_grid_xt<XT>::neighbor_distances_impl(size_type row,
                                                            size_type col) const noexcept
        -> const distance_list&
    {
        return m_neighbor_distances_queen.at(gcode(row, col));
    }

    template <class XT>
    template <raster_connect RC, std::enable_if_t<RC == raster_connect::rook, int>>
    inline auto raster_grid_xt<XT>::neighbor_distances_impl(size_type row,
                                                            size_type col) const noexcept
        -> const distance_list&
    {
        return m_neighbor_distances_rook.at(gcode(row, col));
    }

    template <class XT>
    template <raster_connect RC>
    inline auto raster_grid_xt<XT>::neighbor_distances(size_type row,
                                                    size_type col) const noexcept
        -> const distance_list&
    {
        return neighbor_distances_impl<RC>(row, col);
    }

    template <class XT>
    template <raster_connect RC>
    inline auto raster_grid_xt<XT>::neighbor_indices(size_type row,
                                                    size_type col) const noexcept
        -> neighbor_indices_type
    {
        xt::xtensor_fixed<size_type, xt::xshape<2>> idx {row, col};

        const auto& offsets = neighbor_offsets<RC>(gcode(row, col));
        auto indices = neighbor_indices_type(offsets.shape());

        auto id_it = indices.begin();
        for(auto it=offsets.cbegin(); it!=offsets.cend(); ++it)
        {
            xt::noalias(*id_it++) = *it + idx;
        }

        return indices;
    }

    template <class XT>
    template <raster_connect RC, class E>
    inline auto raster_grid_xt<XT>::neighbor_view(E&& field,
                                                size_type row,
                                                size_type col) const noexcept
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

}