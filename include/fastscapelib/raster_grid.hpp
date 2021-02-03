/**
 * A raster grid represents a 2D uniform grid.
 */
#ifndef FASTSCAPELIB_RASTER_GRID_H
#define FASTSCAPELIB_RASTER_GRID_H

#include "fastscapelib/structured_grid.hpp"
#include "fastscapelib/profile_grid.hpp"


namespace fastscapelib
{
    /*****************
     * Grid elements *
     *****************/

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

        bool is_vertical_looped() const;

    protected:

        void check_looped_symmetrical() const;        
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
    inline bool raster_boundary_status::is_vertical_looped() const
    {
        return is_looped(top) && is_looped(bottom);
    }

    /**
     * Throw an error if the boundary status is not symmetrical.
     */
    inline void raster_boundary_status::check_looped_symmetrical() const
    {
        if (is_looped(left) ^ is_looped(right) || is_looped(top) ^ is_looped(bottom))
        {
            throw std::invalid_argument("looped boundaries are not symmetrical");
        }
    }

    /*******************
    * Raster grid (2D) *
    ********************/

    /**
     * Raster grid node connectivity.
     */
    enum class raster_connect
    {
        rook = 0,       /**< 4-connectivity (vertical or horizontal) */
        queen = 1,      /**< 8-connectivity (including diagonals) */
        bishop = 2      /**< 4-connectivity (only diagonals) */
    };

    /*************************
     * raster_neighbors_base *
     *************************/

    struct raster_neighbors_base
    {
        using neighbors_offsets_type = xt::xtensor<xt::xtensor_fixed<std::ptrdiff_t, xt::xshape<2>>, 1>;
        using neighbors_count_type = std::uint8_t;
    };

    /********************
     * raster_neighbors *
     ********************/

    template <raster_connect RC>
    struct raster_neighbors
    {};


    template <>
    struct raster_neighbors<raster_connect::queen>: public raster_neighbors_base
    {
        // TODO: understand why using a static constexpr raises a linker error on CI
        static constexpr std::uint8_t max_neighbors()
        {
            return 8u;
        }

        neighbors_offsets_type node_neighbors_offsets(std::ptrdiff_t up,
                                                      std::ptrdiff_t down,
                                                      std::ptrdiff_t left,
                                                      std::ptrdiff_t right) const;

        std::array<neighbors_count_type, 9> build_neighbors_count(const raster_boundary_status& status_at_bounds) const;
    };

    /**
     * Given the possible moves (up/down/left/right) depending on grid periodicity,
     * return a list of offsets corresponding to "queen" neighbors relative locations.
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
    inline auto raster_neighbors<raster_connect::queen>::node_neighbors_offsets(std::ptrdiff_t up,
                                                                                std::ptrdiff_t down,
                                                                                std::ptrdiff_t left,
                                                                                std::ptrdiff_t right)
        const -> neighbors_offsets_type
    {
        xt::xtensor<bool, 1> mask {
            (up != 0 && left != 0), up != 0, (up != 0 && right != 0),
            left != 0, right != 0,
            (down != 0 && left != 0), down != 0, (down != 0 && right != 0)
        };

        neighbors_offsets_type full_offsets {
            {up, left}, {up, 0}, {up, right},
            {0, left}, {0, right},
            {down, left}, {down, 0}, {down, right}
        };

        return xt::filter(full_offsets, mask);
    }

    /**
     * Build the neighbors count of each node code (e.g. position on the grid) 
     * depending on the periodicity of the grid.
     */
    inline auto raster_neighbors<raster_connect::queen>::build_neighbors_count(const raster_boundary_status& status_at_bounds)
        const -> std::array<neighbors_count_type, 9>
    {
        std::array<neighbors_count_type, 9> neighbors_count;

        if (status_at_bounds.is_vertical_looped() && status_at_bounds.is_horizontal_looped())
        {
            neighbors_count.fill(8);
        } else 
        if (status_at_bounds.is_vertical_looped())
        {
            neighbors_count = std::array<std::uint8_t, 9>({5, 8, 5, 5, 8, 5, 5, 8, 5});
        } else
        if (status_at_bounds.is_horizontal_looped())
        {
            neighbors_count = std::array<std::uint8_t, 9>({5, 5, 5, 8, 8, 8, 5, 5, 5});
        } else
        {
            neighbors_count = std::array<std::uint8_t, 9>({3, 5, 3, 5, 8, 5, 3, 5, 3});
        }

        return neighbors_count;
    }


    template <>
    struct raster_neighbors<raster_connect::rook>: public raster_neighbors_base
    {
        static constexpr std::uint8_t max_neighbors()
        {
            return 4u;
        }

        neighbors_offsets_type node_neighbors_offsets(std::ptrdiff_t up,
                                                      std::ptrdiff_t down,
                                                      std::ptrdiff_t left,
                                                      std::ptrdiff_t right) const;
        
        std::array<neighbors_count_type, 9> build_neighbors_count(const raster_boundary_status& status_at_bounds) const;
    };

    /**
     * Given the possible moves (up/down/left/right) depending on grid periodicity,
     * return a list of offsets corresponding to "rook" neighbors relative locations.
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
    inline auto raster_neighbors<raster_connect::rook>::node_neighbors_offsets(std::ptrdiff_t up,
                                                                               std::ptrdiff_t down,
                                                                               std::ptrdiff_t left,
                                                                               std::ptrdiff_t right)
        const -> neighbors_offsets_type
    {
        xt::xtensor<bool, 1> mask {
            up != 0, left != 0, right != 0, down != 0
        };

        neighbors_offsets_type full_offsets {
            {up, 0}, {0, left}, {0, right}, {down, 0}
        };

        return xt::filter(full_offsets, mask);
    }

    /**
     * Build the neighbors count of each node code (e.g. position on the grid) 
     * depending on the periodicity of the grid.
     */
    inline auto raster_neighbors<raster_connect::rook>::build_neighbors_count(const raster_boundary_status& status_at_bounds)
        const -> std::array<neighbors_count_type, 9>
    {
        std::array<neighbors_count_type, 9> neighbors_count;

        if (status_at_bounds.is_vertical_looped() && status_at_bounds.is_horizontal_looped())
        {
            neighbors_count.fill(4);
        } else 
        if (status_at_bounds.is_vertical_looped())
        {
            neighbors_count = std::array<std::uint8_t, 9>({3, 4, 3, 3, 4, 3, 3, 4, 3});
        } else
        if (status_at_bounds.is_horizontal_looped())
        {
            neighbors_count = std::array<std::uint8_t, 9>({3, 3, 3, 4, 4, 4, 3, 3, 3});
        } else
        {
            neighbors_count = std::array<std::uint8_t, 9>({2, 3, 2, 3, 4, 3, 2, 3, 2});
        }

        return neighbors_count;
    }


    template <>
    struct raster_neighbors<raster_connect::bishop>: public raster_neighbors_base
    {
        static constexpr std::uint8_t max_neighbors()
        {
            return 4u;
        }

        neighbors_offsets_type node_neighbors_offsets(std::ptrdiff_t up,
                                                      std::ptrdiff_t down,
                                                      std::ptrdiff_t left,
                                                      std::ptrdiff_t right) const;

        std::array<neighbors_count_type, 9> build_neighbors_count(const raster_boundary_status& status_at_bounds) const;
    };

    /**
     * Given the possible moves (up/down/left/right) depending on grid periodicity,
     * return a list of offsets corresponding to "bishop" neighbors relative locations.
     * Offsets are returned in row-major layout when accessible (else dropped).
     *
     * Examples:
     * 
     *         All          w/o left       w/o top
     *      0       1            0
     *        \   /            /
     *          N            N                N
     *        /   \            \            /   \
     *      2       3            1        0       1
     */
    inline auto raster_neighbors<raster_connect::bishop>::node_neighbors_offsets(std::ptrdiff_t up,
                                                                                 std::ptrdiff_t down,
                                                                                 std::ptrdiff_t left,
                                                                                 std::ptrdiff_t right)
        const -> neighbors_offsets_type
    {
        xt::xtensor<bool, 1> mask {
            (up != 0 && left != 0), (up != 0 && right != 0),
            (down != 0 && left != 0), (down != 0 && right != 0)
        };

        neighbors_offsets_type full_offsets {
            {up, left}, {up, right},
            {down, left}, {down, right}
        };

        return xt::filter(full_offsets, mask);
    }

    /**
     * Build the neighbors count of each node code (e.g. position on the grid) 
     * depending on the periodicity of the grid.
     */
    inline auto raster_neighbors<raster_connect::bishop>::build_neighbors_count(const raster_boundary_status& status_at_bounds)
        const -> std::array<neighbors_count_type, 9>
    {
        std::array<neighbors_count_type, 9> neighbors_count;

        if (status_at_bounds.is_vertical_looped() && status_at_bounds.is_horizontal_looped())
        {
            neighbors_count.fill(4);
        } else 
        if (status_at_bounds.is_vertical_looped())
        {
            neighbors_count = std::array<std::uint8_t, 9>({2, 4, 2, 2, 4, 2, 2, 4, 2});
        } else
        if (status_at_bounds.is_horizontal_looped())
        {
            neighbors_count = std::array<std::uint8_t, 9>({2, 2, 2, 4, 4, 4, 2, 2, 2});
        } else
        {
            neighbors_count = std::array<std::uint8_t, 9>({1, 2, 1, 2, 4, 2, 1, 2, 1});
        }

        return neighbors_count;
    }

    /********************
     * grid_inner_types *
     ********************/

    template <class XT, raster_connect RC, class C>
    class raster_grid_xt;
    
    template <class XT, raster_connect RC, class C>
    struct grid_inner_types<raster_grid_xt<XT, RC, C>>
    {
        static constexpr bool is_structured = true;
        static constexpr bool is_uniform = true;
        static constexpr uint8_t max_neighbors = raster_neighbors<RC>::max_neighbors();
        static constexpr std::size_t xt_ndims = 2;
        
        using xt_selector = XT;
        using xt_type = xt_container_t<xt_selector, int, xt_ndims>;

        using size_type = typename xt_type::size_type;
        using shape_type = typename xt_type::shape_type;
        using distance_type = double;
        using length_type = xt::xtensor_fixed<distance_type, xt::xshape<2>>;
        using spacing_type = length_type;

        using code_type = std::uint8_t;

        using neighbors_count_type = typename raster_neighbors_base::neighbors_count_type;
        using neighbors_distances_impl_type = typename std::array<distance_type, max_neighbors>;

        using boundary_status_type = raster_boundary_status;
        using node_status_type = xt_container_t<xt_selector, node_status, xt_ndims>;        
    };


    /**
     * @class raster_grid_xt
     * @brief 2-dimensional uniform (raster) grid.
     *
     * @tparam XT xtensor container selector for data array members.
     */
    template <class XT, raster_connect RC, class C = detail::neighbors_cache<raster_neighbors<RC>::max_neighbors()>>
    class raster_grid_xt : public structured_grid<raster_grid_xt<XT, RC, C>, C>,
                           public raster_neighbors<RC>
    {
    public:

        using self_type = raster_grid_xt<XT, RC, C>;
        using base_type = structured_grid<self_type, C>;
        using inner_types = grid_inner_types<self_type>;

        using xt_selector = typename inner_types::xt_selector;
        static constexpr std::size_t xt_ndims = inner_types::xt_ndims;
        static constexpr std::uint8_t max_neighbors(){ return inner_types::max_neighbors; };

        using size_type = typename inner_types::size_type;
        using shape_type = typename inner_types::shape_type;
        using length_type = typename inner_types::length_type;
        using spacing_type = typename inner_types::spacing_type;
        using distance_type = typename inner_types::distance_type;

        using code_type = typename inner_types::code_type;

        using neighbors_type = typename base_type::neighbors_type;
        using neighbors_count_type = typename inner_types::neighbors_count_type;
        using neighbors_indices_type = typename base_type::neighbors_indices_type;
        using neighbors_distances_type = typename base_type::neighbors_distances_type;
        using neighbors_indices_raster_type = xt::xtensor<xt::xtensor_fixed<size_type, xt::xshape<2>>, 1>;

        using boundary_status_type = typename inner_types::boundary_status_type;
        using node_status_type = typename inner_types::node_status_type;

        raster_grid_xt(const shape_type& shape,
                       const spacing_type& spacing,
                       const boundary_status_type& status_at_bounds,
                       const std::vector<raster_node>& status_at_nodes = {});

        static raster_grid_xt from_length(const shape_type& shape,
                                          const length_type& length,
                                          const boundary_status_type& status_at_bounds,
                                          const std::vector<raster_node>& status_at_nodes = {});

        template<raster_connect RCA, class CA = C>
        static raster_grid_xt<XT, RCA, CA> clone(const raster_grid_xt& grid);

        using base_type::neighbors_distances;

        using base_type::neighbors_indices;
        inline neighbors_indices_raster_type neighbors_indices(const size_type& row,
                                                               const size_type& col) const noexcept;

        template <class E>
        inline auto neighbor_view(E&& field, size_type idx) const noexcept;

        code_type node_code(const size_type& row, const size_type& col) const noexcept;
        code_type node_code(const size_type& idx) const noexcept;

        inline const neighbors_count_type& neighbors_count(const size_type& idx) const noexcept;
        inline const neighbors_count_type& neighbors_count(const code_type& code) const noexcept;

        shape_type shape() const noexcept;

        boundary_status_type status_at_bounds() const noexcept;

    private:

        using neighbors_distances_impl_type = typename inner_types::neighbors_distances_impl_type;
        using neighbors_indices_impl_type = typename base_type::neighbors_indices_impl_type;
        using neighbors_offsets_type = typename raster_neighbors<RC>::neighbors_offsets_type;

        using coded_ncount_type = std::array<neighbors_count_type, 9>;
        using coded_noffsets_type = std::array<neighbors_offsets_type, 9>;        
        using coded_ndistances_type = std::array<neighbors_distances_impl_type, 9>;
        using node_code_type = xt::xtensor<code_type, 1>;

        shape_type m_shape;
        size_type m_size;
        spacing_type m_spacing;
        length_type m_length;

        node_status_type m_status_at_nodes;
        boundary_status_type m_status_at_bounds;

        struct corner_node
        {
            size_type row;
            size_type col;
            node_status row_border;
            node_status col_border;
        };

        node_code_type m_nodes_codes;

        coded_ncount_type m_neighbors_count;
        coded_noffsets_type m_neighbor_offsets;
        coded_ndistances_type m_neighbor_distances;

        void set_status_at_nodes(const std::vector<raster_node>& status_at_nodes);

        void build_nodes_codes();
        coded_noffsets_type build_coded_neighbors_offsets();
        coded_ndistances_type build_coded_neighbors_distances();

        inline const neighbors_offsets_type& neighbor_offsets(code_type code) const noexcept;

        inline neighbors_indices_impl_type get_neighbors_indices(const size_type& idx) const noexcept;

        inline const neighbors_distances_impl_type& neighbors_distances_impl(const size_type& idx) const noexcept;

        void neighbors_indices_impl(neighbors_indices_impl_type& neighbors, const size_type& idx) const;

        friend class structured_grid<self_type, C>;
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
    template <class XT, raster_connect RC, class C>
    raster_grid_xt<XT, RC, C>::raster_grid_xt(const shape_type& shape,
                                              const spacing_type& spacing,
                                              const boundary_status_type& status_at_bounds,
                                              const std::vector<raster_node>& status_at_nodes)
        : base_type(shape[0] * shape[1]), m_shape(shape), m_spacing(spacing), m_status_at_bounds(status_at_bounds)
    {
        m_size = shape[0] * shape[1];
        m_length = (xt::adapt(shape) - 1) * spacing;

        build_nodes_codes();
        m_neighbors_count = this->build_neighbors_count(status_at_bounds);        

        set_status_at_nodes(status_at_nodes);

        m_neighbor_offsets = build_coded_neighbors_offsets();
        m_neighbor_distances = build_coded_neighbors_distances();
    }
    //@}

    /**
     * @name Factories
     */
    //@{
   /**
     * Creates a new raster grid.
     *
     * @param shape Shape of the grid (number of rows and cols).
     * @param length Total physical lengths of the grid.
     * @param status_at_bounds Status at boundary nodes (left & right grid edges).
     * @param status_at_nodes Manually define the status at any node on the grid.
     */
    template <class XT, raster_connect RC, class C>
    raster_grid_xt<XT, RC, C> raster_grid_xt<XT, RC, C>::from_length(const shape_type& shape,
                                                             const length_type& length,
                                                             const boundary_status_type& status_at_bounds,
                                                             const std::vector<raster_node>& status_at_nodes)
    {
        spacing_type spacing = length / (xt::adapt(shape) - 1);
        return raster_grid_xt<XT, RC, C>(shape, spacing, status_at_bounds, status_at_nodes);
    }
    //@}

    template <class XT, raster_connect RC, class C>
    template<raster_connect RCA, class CA>
    raster_grid_xt<XT, RCA, CA> raster_grid_xt<XT, RC, C>::clone(const raster_grid_xt& grid)
    {
        return raster_grid_xt<XT, RCA, CA>(grid.shape(), grid.spacing(), grid.status_at_bounds());
    }

    template <class XT, raster_connect RC, class C>
    auto raster_grid_xt<XT, RC, C>::shape() const noexcept
        -> shape_type
    {
        return m_shape;
    }

    template <class XT, raster_connect RC, class C>
    auto raster_grid_xt<XT, RC, C>::status_at_bounds() const noexcept
        -> boundary_status_type
    {
        return m_status_at_bounds;
    }

    template <class XT, raster_connect RC, class C>
    void raster_grid_xt<XT, RC, C>::set_status_at_nodes(const std::vector<raster_node>& status_at_nodes)
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
    template <class XT, raster_connect RC, class C>
    void raster_grid_xt<XT, RC, C>::build_nodes_codes()
    {
        std::array<std::vector<code_type>, 2> gcode_rc;

        for (std::uint8_t dim=0; dim<2; ++dim)
        {
            auto fill_value = static_cast<std::uint8_t>(3 - dim * 2);
            std::vector<std::uint8_t> gcode_component(m_shape[dim], fill_value);

            gcode_component[0] = 0;
            gcode_component[m_shape[dim]-1] = static_cast<std::uint8_t>(fill_value * 2);

            gcode_rc[dim] = gcode_component;
        }

        m_nodes_codes.resize({m_size});
        for (std::size_t r=0; r<m_shape[0]; ++r)
        {
            for (std::size_t c=0; c<m_shape[1]; ++c)
            {
                m_nodes_codes[r*m_shape[1]+c] = static_cast<std::uint8_t>(gcode_rc[0][r] + gcode_rc[1][c]);
            }
        }
    }

    template <class XT, raster_connect RC, class C>
    auto raster_grid_xt<XT, RC, C>::neighbors_count(const size_type& idx) const noexcept
    -> const neighbors_count_type&
    {
        return m_neighbors_count[m_nodes_codes(idx)];
    }

    template <class XT, raster_connect RC, class C>
    auto raster_grid_xt<XT, RC, C>::neighbors_count(const code_type& code) const noexcept
    -> const neighbors_count_type&
    {
        return m_neighbors_count[code];
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
    template <class XT, raster_connect RC, class C>
    inline auto raster_grid_xt<XT, RC, C>::node_code(const size_type& row, const size_type& col) const noexcept
    -> code_type
    {
        return m_nodes_codes[row * m_shape[1] + col];
    }

    template <class XT, raster_connect RC, class C>
    inline auto raster_grid_xt<XT, RC, C>::node_code(const size_type& idx) const noexcept
    -> code_type
    {
        return m_nodes_codes[idx];
    }

    /**
     * Pre-store the (row, col) index offsets of grid node neighbors for
     * each of the 9-characteristic locations on the grid.
     *
     * Those offsets take into account looped boundary conditions (if any).
     * The order of the returned offsets corresponds to the row-major layout.
     */
    template <class XT, raster_connect RC, class C>
    auto raster_grid_xt<XT, RC, C>::build_coded_neighbors_offsets()
        -> coded_noffsets_type
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
            this->node_neighbors_offsets(dr, 1, dc, 1),      // top-left corner
            this->node_neighbors_offsets(dr, 1, -1, 1),      // top border
            this->node_neighbors_offsets(dr, 1, -1, -dc),    // top-right corner
            this->node_neighbors_offsets(-1, 1, dc, 1),      // left border
            this->node_neighbors_offsets(-1, 1, -1, 1),      // inner
            this->node_neighbors_offsets(-1, 1, -1, -dc),    // right border
            this->node_neighbors_offsets(-1, -dr, dc, 1),    // bottom-left corner
            this->node_neighbors_offsets(-1, -dr, -1, 1),    // bottom border
            this->node_neighbors_offsets(-1, -dr, -1, -dc),  // bottom-right corner
        };
    }

    template <class XT, raster_connect RC, class C>
    inline auto raster_grid_xt<XT, RC, C>::neighbor_offsets(code_type code) const noexcept
        -> const neighbors_offsets_type&
    {
        return m_neighbor_offsets[code];
    }

    template <class XT, raster_connect RC, class C>
    auto raster_grid_xt<XT, RC, C>::build_coded_neighbors_distances()
        -> coded_ndistances_type
    {
        coded_ndistances_type nb_distances;
        auto xspacing = m_spacing;

        auto to_dist = [&xspacing](auto&& offset) -> double {
                        auto drc = xt::where(xt::equal(offset, 0), 0., 1.) * xspacing;
                        return std::sqrt(xt::sum(xt::square(drc))(0));
                    };

        for(std::uint8_t k=0; k<9; ++k)
        {
            auto offsets = neighbor_offsets(k);
            auto distances = neighbors_distances_impl_type();

            std::transform(offsets.cbegin(), offsets.cend(), distances.begin(), to_dist);

            nb_distances[k] = distances;
        }

        return nb_distances;
    }

    template <class XT, raster_connect RC, class C>
    inline auto raster_grid_xt<XT, RC, C>::neighbors_distances_impl(const size_type& idx) const noexcept
        -> const neighbors_distances_impl_type&
    {
        return m_neighbor_distances[node_code(idx)];
    }

    template <class XT, raster_connect RC, class C>
    inline auto raster_grid_xt<XT, RC, C>::neighbors_indices(const size_type& row,
                                                             const size_type& col) const noexcept
        -> neighbors_indices_raster_type
    {
        xt::xtensor_fixed<size_type, xt::xshape<2>> idx {row, col};

        const auto& offsets = neighbor_offsets(node_code(row, col));
        auto indices = neighbors_indices_raster_type(offsets.shape());

        auto id_it = indices.begin();
        for(auto it=offsets.cbegin(); it!=offsets.cend(); ++it)
        {
            xt::noalias(*id_it++) = *it + idx;
        }

        return indices;
    }

    template <class XT, raster_connect RC, class C>
    inline auto raster_grid_xt<XT, RC, C>::get_neighbors_indices(const size_type& idx) const noexcept
        -> neighbors_indices_impl_type
    {
        const auto& offsets = neighbor_offsets(node_code(idx));
        neighbors_indices_impl_type indices;

        auto id_it = indices.begin();
        for(auto it=offsets.cbegin(); it!=offsets.cend(); ++it)
        {
            (*id_it++) = static_cast<size_type>((*it)[0]) * m_shape[1] + static_cast<size_type>((*it)[1]) + idx;
        }

        return indices;
    }

    template <class XT, raster_connect RC, class C>    
    inline auto raster_grid_xt<XT, RC, C>::neighbors_indices_impl(neighbors_indices_impl_type& neighbors, const size_type& idx) const
        -> void
    {
        auto indices = get_neighbors_indices(idx);
        
        for (size_type i=0; i<indices.size(); ++i)
        {
            neighbors.at(i) = indices[i];
        }
    }

    template <class XT, raster_connect RC, class C>
    template <class E>
    inline auto raster_grid_xt<XT, RC, C>::neighbor_view(E&& field,
                                                  size_type idx) const noexcept
    {
        return xt::index_view(std::forward<E>(field), get_neighbors_indices(idx));
    }

    /**
     * @typedef raster_grid
     * Alias template on raster_grid_xt with ``xt::xtensor`` used
     * as array container type for data members.
     *
     * This is mainly for convenience when using in C++ applications.
     *
     */
    using raster_grid = raster_grid_xt<xtensor_selector, raster_connect::queen>;

}

#endif