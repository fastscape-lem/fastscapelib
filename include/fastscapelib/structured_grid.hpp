/**
 * grids (channel, raster, etc.) hold and manage grid/mesh geometry,
 * topology and boundary conditions.
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <type_traits>
#include <map>
#include <vector>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xfixed.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xindex_view.hpp"
#include "xtensor/xnoalias.hpp"

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
    enum class node_status: std::uint8_t
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
    protected:

        bool is_looped(node_status status);
    };

    inline bool boundary_status::is_looped(node_status status)
    {
        return status == node_status::looped_boundary;
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
     * Neighbor node.
     */
    struct neighbor
    {
        std::size_t idx;     /**< Index of the neighbor node */
        double distance;     /**< Distance to the neighbor node */
        node_status status;  /**< Status at the neighbor node */

        bool operator==(const neighbor& rhs) const
        {
            return (idx == rhs.idx) && (distance == rhs.distance) && (status == rhs.status);
        }
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

    //***************************
    //* Structured grid interface
    //***************************

    template <class G>
    struct grid_inner_types;

    /**
     * @class structured_grid
     * @brief Add a common interface to all structured grid types.
     *
     * This class only defines a basic interface for all structured grid types.
     * It does not embed any data member, this responsibility
     * is delegated to the inheriting classes.
     *
     * @tparam G Derived grid type.
     */
    template <class G>
    class structured_grid
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
        structured_grid() = default;
        ~structured_grid() = default;

        const derived_grid_type& derived_grid() const & noexcept;
    };

    template <class G>
    inline auto structured_grid<G>::derived_grid() const & noexcept -> const derived_grid_type&
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
    inline auto structured_grid<G>::size() const noexcept -> size_type
    {
        return derived_grid().m_size;
    }

    /**
     * Returns the (uniform) spacing between two adjacent grid nodes.
     *
     * Returns a copy of the value for 1-d grids or a constant reference otherwise.
     */
    template <class G>

    inline auto structured_grid<G>::spacing() const noexcept -> spacing_t
    {
        return derived_grid().m_spacing;
    }

    /**
     * Returns a constant reference to the array of status at grid nodes.
     */
    template <class G>
    inline auto structured_grid<G>::status_at_nodes() const -> const node_status_type&
    {
        return derived_grid().m_status_at_nodes;
    }

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
    inline auto structured_grid<G>::neighbors(std::size_t idx) const -> const neighbors_type&
    {
        return derived_grid().neighbors_impl(idx);
    }
    //@}
}