/**
 * A profile grid is a one dimensional grid.
 */
#pragma once

#include "fastscapelib/structured_grid.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <map>
#include <vector>


namespace fastscapelib
{

    /**
     * Status at grid boundary nodes.
     */
    class profile_boundary_status: public boundary_status
    {
    public:

        node_status left = node_status::core;    /**< Status at left edge/border node(s) */
        node_status right = node_status::core;   /**< Status at right edge/border node(s) */

        profile_boundary_status(const node_status status);
        profile_boundary_status(const node_status left_status, const node_status right_status);
        profile_boundary_status(const std::array<node_status, 2>& status);

        bool is_horizontal_looped();

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
        : left(status), right(status)
    {
        check_looped_symmetrical();
    }

    /**
     * Set the same status for all boundary nodes.
     */
    inline profile_boundary_status::profile_boundary_status(node_status left_status, node_status right_status)
        : left(left_status), right(right_status)
    {
        check_looped_symmetrical();
    }

    /**
     * Set status at the left and right edge nodes of a profile grid.
     */
    inline profile_boundary_status::profile_boundary_status(const std::array<node_status, 2>& status)
        : left(status[0]), right(status[1])
    {
        check_looped_symmetrical();
    }
    //@}

    /**
     * Return true if periodic (looped) conditions are set for ``left`` and ``right``.
     */
    inline bool profile_boundary_status::is_horizontal_looped()
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

        using boundary_status_type = profile_boundary_status;
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
    class profile_grid_xt : public structured_grid<profile_grid_xt<XT>>
    {
    public:

        using self_type = profile_grid_xt<XT>;
        using base_type = structured_grid<self_type>;
        using inner_types = grid_inner_types<self_type>;

        using xt_selector = typename inner_types::xt_selector;
        using xt_ndims = typename inner_types::xt_ndims;

        using size_type = typename inner_types::size_type;
        using shape_type = typename inner_types::shape_type;

        using spacing_type = typename inner_types::spacing_type;
        using neighbors_type = typename inner_types::neighbors_type;

        using boundary_status_type = typename inner_types::boundary_status_type;
        using node_status_type = typename inner_types::node_status_type;

        profile_grid_xt(size_type size,
                        spacing_type spacing,
                        const boundary_status_type& status_at_bounds,
                        const std::vector<node>& status_at_nodes = {});

    protected:

        shape_type m_shape;
        size_type m_size;
        spacing_type m_spacing;

        node_status_type m_status_at_nodes;
        boundary_status_type m_status_at_bounds;
        void set_status_at_nodes(const std::vector<node>& status_at_nodes);

        static constexpr std::array<std::ptrdiff_t, 3> offsets { {0, -1, 1} };
        std::vector<neighbors_type> m_all_neighbors;
        void precompute_neighbors();

        const neighbors_type& neighbors_impl(std::size_t idx) const;

        friend class structured_grid<self_type>;
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
                                        const boundary_status_type& status_at_bounds,
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
            else if (temp_status_at_nodes.at(n.idx) == node_status::looped_boundary)
            {
                throw std::invalid_argument("cannot overwrite the status of a "
                                            "looped boundary node");
            }

            temp_status_at_nodes.at(n.idx) = n.status;
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
        return m_all_neighbors.at(idx);
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