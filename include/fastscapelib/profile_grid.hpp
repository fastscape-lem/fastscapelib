/**
 * A profile grid is a one dimensional grid.
 */
#ifndef FASTSCAPELIB_PROFILE_GRID_H
#define FASTSCAPELIB_PROFILE_GRID_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <map>
#include <vector>

#include "fastscapelib/grid.hpp"
#include "fastscapelib/structured_grid.hpp"


namespace fastscapelib
{

    /**
     * Status at grid boundary nodes.
     */
    class profile_boundary_status : public boundary_status
    {
    public:
        node_status left = node_status::core;  /**< Status at left edge/border node(s) */
        node_status right = node_status::core; /**< Status at right edge/border node(s) */

        profile_boundary_status(const node_status status);
        profile_boundary_status(const node_status left_status, const node_status right_status);
        profile_boundary_status(const std::array<node_status, 2>& status);

        bool is_horizontal_looped() const;

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
    inline profile_boundary_status::profile_boundary_status(node_status status)
        : left(status)
        , right(status)
    {
        check_looped_symmetrical();
    }

    /**
     * Set status at the left and right edge nodes of a profile grid.
     */
    inline profile_boundary_status::profile_boundary_status(node_status left_status,
                                                            node_status right_status)
        : left(left_status)
        , right(right_status)
    {
        check_looped_symmetrical();
    }

    /**
     * Set status at the left and right edge nodes of a profile grid.
     */
    inline profile_boundary_status::profile_boundary_status(
        const std::array<node_status, 2>& status)
        : left(status[0])
        , right(status[1])
    {
        check_looped_symmetrical();
    }
    //@}

    /**
     * Return true if periodic (looped) conditions are set for ``left`` and ``right``.
     */
    inline bool profile_boundary_status::is_horizontal_looped() const
    {
        return is_looped(left) && is_looped(right);
    }

    inline void profile_boundary_status::check_looped_symmetrical()
    {
        if (is_looped(left) ^ is_looped(right))
        {
            throw std::invalid_argument("looped boundaries are not symmetrical");
        }
    }

    //*******************
    //* Profile grid (1D)
    //*******************

    template <class XT, class C>
    class profile_grid_xt;

    template <class XT, class C>
    struct grid_inner_types<profile_grid_xt<XT, C>>
    {
        static constexpr bool is_structured = true;
        static constexpr bool is_uniform = true;
        static constexpr uint8_t max_neighbors = 2;
        static constexpr std::size_t xt_ndims = 1;

        using xt_selector = XT;
        using xt_type = xt_container_t<xt_selector, int, xt_ndims>;

        using size_type = typename xt_type::size_type;
        using shape_type = typename xt_type::shape_type;
        using length_type = double;
        using distance_type = double;
        using spacing_type = length_type;

        using code_type = std::uint8_t;

        using neighbors_count_type = std::uint8_t;
        using neighbors_distances_impl_type = typename std::array<distance_type, max_neighbors>;

        using boundary_status_type = profile_boundary_status;
        using node_status_type = xt_container_t<xt_selector, node_status, xt_ndims>;
    };

    /**
     * @class profile_grid_xt
     * @brief 1-dimensional uniform grid.
     *
     * Used for modeling single channel or hillslope profiles.
     *
     * @tparam XT xtensor container selector for data array members.
     * @tparam C Grid neighbors cache type.
     */
    template <class XT, class C = neighbors_cache<2>>
    class profile_grid_xt : public structured_grid<profile_grid_xt<XT, C>, C>
    {
    public:
        using self_type = profile_grid_xt<XT, C>;
        using base_type = structured_grid<self_type, C>;
        using neighbors_cache_type = C;
        using inner_types = grid_inner_types<self_type>;

        using xt_selector = typename inner_types::xt_selector;
        static constexpr std::size_t xt_ndims = inner_types::xt_ndims;
        static constexpr std::uint8_t max_neighbors()
        {
            return inner_types::max_neighbors;
        };

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

        using boundary_status_type = typename inner_types::boundary_status_type;
        using node_status_type = typename inner_types::node_status_type;

        profile_grid_xt(size_type size,
                        spacing_type spacing,
                        const boundary_status_type& status_at_bounds,
                        const std::vector<node>& status_at_nodes = {});

        static profile_grid_xt from_length(size_type size,
                                           length_type length,
                                           const boundary_status_type& status_at_bounds,
                                           const std::vector<node>& status_at_nodes = {});

    protected:
        using neighbors_distances_impl_type = typename inner_types::neighbors_distances_impl_type;
        using neighbors_indices_impl_type = typename base_type::neighbors_indices_impl_type;

        using coded_ndistances_type = std::array<neighbors_distances_impl_type, 3>;

        shape_type m_shape;
        size_type m_size;
        spacing_type m_spacing;
        length_type m_length;
        distance_type m_node_area;

        xt::xtensor<code_type, 1> m_gcode_idx;

        node_status_type m_status_at_nodes;
        boundary_status_type m_status_at_bounds;

        void set_status_at_nodes(const std::vector<node>& status_at_nodes);

        static constexpr std::array<std::ptrdiff_t, 3> offsets{ { 0, -1, 1 } };
        std::vector<neighbors_type> m_all_neighbors;

        coded_ndistances_type m_neighbors_distances;
        void build_neighbors_distances();

        void build_gcode();
        code_type gcode(const size_type& idx) const;

        std::array<neighbors_count_type, 3> m_neighbors_count;
        void build_neighbors_count();
        inline const neighbors_count_type& neighbors_count(const size_type& idx) const noexcept;
        inline const neighbors_count_type& neighbors_count(const code_type& code) const noexcept;

        void neighbors_indices_impl(neighbors_indices_impl_type& neighbors,
                                    const size_type& idx) const;

        const neighbors_distances_impl_type& neighbors_distances_impl(const size_type& idx) const;

        friend class structured_grid<self_type, neighbors_cache_type>;
        friend class grid<self_type, neighbors_cache_type>;
    };

    template <class XT, class C>
    constexpr std::array<std::ptrdiff_t, 3> profile_grid_xt<XT, C>::offsets;

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
    template <class XT, class C>
    profile_grid_xt<XT, C>::profile_grid_xt(size_type size,
                                            spacing_type spacing,
                                            const boundary_status_type& status_at_bounds,
                                            const std::vector<node>& status_at_nodes)
        : base_type(size)
        , m_size(size)
        , m_spacing(spacing)
        , m_status_at_bounds(status_at_bounds)
    {
        m_shape = { static_cast<typename shape_type::value_type>(m_size) };
        m_length = static_cast<spacing_type>(size - 1) * spacing;
        m_node_area = spacing;
        set_status_at_nodes(status_at_nodes);
        build_gcode();
        build_neighbors_distances();
        build_neighbors_count();
    }
    //@}

    /**
     * @name Factories
     */
    //@{
    /**
     * Creates a new profile grid.
     *
     * @param size Total number of grid nodes.
     * @param length Total physical length of the grid.
     * @param status_at_bounds Status at boundary nodes (left & right grid edges).
     * @param status_at_nodes Manually define the status at any node on the grid.
     */
    template <class XT, class C>
    profile_grid_xt<XT, C> profile_grid_xt<XT, C>::from_length(
        size_type size,
        length_type length,
        const boundary_status_type& status_at_bounds,
        const std::vector<node>& status_at_nodes)
    {
        spacing_type spacing = length / static_cast<length_type>(size - 1);
        return profile_grid_xt<XT, C>(size, spacing, status_at_bounds, status_at_nodes);
    }
    //@}

    template <class XT, class C>
    void profile_grid_xt<XT, C>::build_gcode()
    {
        m_gcode_idx.resize({ m_size });

        m_gcode_idx.fill(1);
        m_gcode_idx[0] = 0;
        m_gcode_idx[m_size - 1] = 2;
    }

    template <class XT, class C>
    auto profile_grid_xt<XT, C>::gcode(const size_type& idx) const -> code_type
    {
        return m_gcode_idx[idx];
    }

    template <class XT, class C>
    void profile_grid_xt<XT, C>::build_neighbors_distances()
    {
        m_neighbors_distances.fill({ m_spacing, m_spacing });
    }

    template <class XT, class C>
    void profile_grid_xt<XT, C>::build_neighbors_count()
    {
        if (m_status_at_bounds.is_horizontal_looped())
        {
            m_neighbors_count = std::array<neighbors_count_type, 3>({ 2, 2, 2 });
        }
        else
        {
            m_neighbors_count = std::array<neighbors_count_type, 3>({ 1, 2, 1 });
        }
    }

    template <class XT, class C>
    auto profile_grid_xt<XT, C>::neighbors_count(const size_type& idx) const noexcept
        -> const neighbors_count_type&
    {
        return m_neighbors_count[gcode(idx)];
    }

    template <class XT, class C>
    auto profile_grid_xt<XT, C>::neighbors_count(const code_type& code) const noexcept
        -> const neighbors_count_type&
    {
        return m_neighbors_count[code];
    }

    template <class XT, class C>
    void profile_grid_xt<XT, C>::set_status_at_nodes(const std::vector<node>& status_at_nodes)
    {
        node_status_type temp_status_at_nodes(m_shape, node_status::core);

        temp_status_at_nodes(0) = m_status_at_bounds.left;
        temp_status_at_nodes(m_size - 1) = m_status_at_bounds.right;

        for (const node& n : status_at_nodes)
        {
            if (n.status == node_status::looped_boundary)
            {
                throw std::invalid_argument("node_status::looped_boundary is not allowed in "
                                            "'status_at_nodes' "
                                            "(use 'status_at_bounds' instead)");
            }
            else if (temp_status_at_nodes.at(n.idx) == node_status::looped_boundary)
            {
                throw std::invalid_argument("cannot overwrite the status of a "
                                            "looped boundary node");
            }

            temp_status_at_nodes.at(n.idx) = n.status;
        }

        m_status_at_nodes = temp_status_at_nodes;
    }

    template <class XT, class C>
    auto profile_grid_xt<XT, C>::neighbors_distances_impl(const size_type& /*idx*/) const
        -> const neighbors_distances_impl_type&
    {
        return m_neighbors_distances[0];
    }

    template <class XT, class C>
    inline auto profile_grid_xt<XT, C>::neighbors_indices_impl(
        neighbors_indices_impl_type& neighbors, const size_type& idx) const -> void
    {
        if (idx == 0)
        {
            if (m_status_at_bounds.is_horizontal_looped())
            {
                neighbors[0] = m_size - 1;
                neighbors[1] = 1;
            }
            else
            {
                neighbors[0] = 1;
            }
        }
        else if (idx == m_size - 1)
        {
            neighbors[0] = m_size - 2;
            if (m_status_at_bounds.is_horizontal_looped())
            {
                neighbors[1] = 0;
            }
        }
        else
        {
            for (std::size_t k = 1; k < 3; ++k)
            {
                std::size_t nb_idx = detail::add_offset(idx, offsets[k]);
                neighbors[k - 1] = nb_idx;
            }
        }
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
}

#endif
