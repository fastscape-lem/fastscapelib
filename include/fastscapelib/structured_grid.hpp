/**
 * grids (channel, raster, etc.) hold and manage grid/mesh geometry,
 * topology and boundary conditions.
 */
#ifndef FASTSCAPELIB_STRUCTURED_GRID_H
#define FASTSCAPELIB_STRUCTURED_GRID_H

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

        bool is_looped(node_status status) const;
    };

    inline bool boundary_status::is_looped(const node_status status) const
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

        /**
         * Provides a cache for neighbors indices queries of a particular node of the grid.
         */
        template <std::uint8_t N>
        class neighbors_cache
        {
        public:
            static constexpr unsigned int cache_width = N;
            using neighbors_indices_type = std::array<std::size_t, N>;

            neighbors_cache(std::size_t size)
                : m_cache(cache_shape_type({size}))
            {
                for (std::size_t i=0; i<size; ++i)
                {
                    m_cache[i].fill(std::numeric_limits<std::size_t>::max());
                }
            }

            bool has(const std::size_t& idx) const
            {
                return m_cache[idx][0] == std::numeric_limits<std::size_t>::max() ? false : true;
            }

            neighbors_indices_type& get(const std::size_t& idx)
            {
                return m_cache[idx];
            }

            neighbors_indices_type& get_storage(const std::size_t& idx)
            {
                return m_cache[idx];
            }

            void store(const std::size_t& idx, const neighbors_indices_type& neighbors_indices)
            {
                m_cache[idx] = neighbors_indices;
            }

            std::size_t cache_size() const
            {
                return m_cache.size();
            };

            std::size_t cache_used() const
            {
                std::size_t count = 0;
                
                for (std::size_t i=0; i < m_cache.size(); ++i)
                {
                    if (m_cache[i][0] != std::numeric_limits<std::size_t>::max()) { count += 1; }
                }

                return count; 
            };

            void reset()
            {
                for (std::size_t i=0; i < m_cache.size(); ++i)
                {
                    m_cache[i].fill(std::numeric_limits<std::size_t>::max());
                }
            }

            void remove(const std::size_t& idx)
            {
                m_cache[idx].fill(std::numeric_limits<std::size_t>::max());
            }

        protected:

            using cache_type = xt::xtensor<neighbors_indices_type, 1>;
            using cache_shape_type = typename cache_type::shape_type;

            cache_type m_cache;
        };
 
 
        /**
         * A pass-through/no cache for neighbors indices queries of a particular node of the grid.
         */
        template <unsigned int N>
        class neighbors_no_cache
        {
        public:
            static constexpr unsigned int cache_width = N;
            using neighbors_indices_type = std::array<std::size_t, N>;

            neighbors_no_cache(std::size_t /*size*/) {}

            bool has(const std::size_t& /*idx*/) const { return false; }

            neighbors_indices_type& get(const std::size_t& /*idx*/)
            {
                return m_node_neighbors;
            }

            neighbors_indices_type& get_storage(const std::size_t& /*idx*/)
            {
                return m_node_neighbors;
            }

            void store(const std::size_t& /*idx*/, const neighbors_indices_type neighbors_indices)
            {
                m_node_neighbors = neighbors_indices;
            }

            std::size_t cache_size() const { return 0; }

            std::size_t cache_used() const { return 0; }

            void reset() {}

            void remove(const std::size_t& /*idx*/) {}

        protected:

            neighbors_indices_type m_node_neighbors;
        };
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
    template <class G, class C>
    class structured_grid
    {
    public:

        using derived_grid_type = G;
        using inner_types = grid_inner_types<derived_grid_type>;
        using neighbors_indices_cache_type = C;

        static_assert(neighbors_indices_cache_type::cache_width >= inner_types::max_neighbors,
                      "Cache width is too small!");

        using size_type = typename inner_types::size_type;
        using shape_type = typename inner_types::shape_type;
        using length_type = typename inner_types::length_type;
        using spacing_type = typename inner_types::spacing_type;
        using distance_type = typename inner_types::distance_type;

        using spacing_t = std::conditional_t<std::is_arithmetic<spacing_type>::value,
                                            spacing_type,
                                            const spacing_type&>;

        using neighbors_type = typename xt::xtensor<neighbor, 1>;
        using neighbors_count_type = typename inner_types::neighbors_count_type;
        using neighbors_indices_type = xt::xtensor<size_type, 1>;
        using neighbors_distances_type = xt::xtensor<distance_type, 1>;

        using node_status_type = typename inner_types::node_status_type;
        
        size_type size() const noexcept;

        spacing_t spacing() const noexcept;

        length_type length() const noexcept;

        const node_status_type& status_at_nodes() const;

        const neighbors_count_type& neighbors_count(const size_type& idx) const;

        neighbors_distances_type neighbors_distances(const size_type& idx) const;

        neighbors_indices_type neighbors_indices(const size_type& idx);

        void neighbors_indices(const size_type& idx, neighbors_indices_type& neighbors_indices);

        neighbors_type neighbors(const size_type& idx);

        void neighbors(const size_type& idx, neighbors_type& neighbors);

        const neighbors_indices_cache_type& neighbors_indices_cache() { return m_neighbors_indices_cache; };

    protected:

        using neighbors_indices_impl_type = typename C::neighbors_indices_type;
        using neighbors_distances_impl_type = typename inner_types::neighbors_distances_impl_type;

        structured_grid(std::size_t size)
            : m_neighbors_indices_cache(neighbors_indices_cache_type(size))
        {};

        ~structured_grid() = default;

        neighbors_indices_cache_type m_neighbors_indices_cache;

        const derived_grid_type& derived_grid() const noexcept;
        derived_grid_type& derived_grid() noexcept;

        const neighbors_indices_impl_type& neighbors_indices_impl(const size_type& idx);

        const neighbors_distances_impl_type& neighbors_distances_impl(const size_type& idx) const;
    };

    template <class G, class C>
    inline auto structured_grid<G, C>::derived_grid() const noexcept
        -> const derived_grid_type&
    {
        return *static_cast<const derived_grid_type*>(this);
    }

    template <class G, class C>
    inline auto structured_grid<G, C>::derived_grid() noexcept
        -> derived_grid_type&
    {
        return *static_cast<derived_grid_type*>(this);
    }

    /**
     * @name Grid properties
     */
    //@{
    /**
     * Returns the total number of grid nodes.
     */
    template <class G, class C>
    inline auto structured_grid<G, C>::size() const noexcept
        -> size_type
    {
        return derived_grid().m_size;
    }

    /**
     * Returns the (uniform) spacing between two adjacent grid nodes.
     *
     * Returns a copy of the value for 1-d grids or a constant reference otherwise.
     */
    template <class G, class C>
    inline auto structured_grid<G, C>::spacing() const noexcept
        -> spacing_t
    {
        return derived_grid().m_spacing;
    }

    /**
     * Returns the length of the grid for all its dimensions.
     */
    template <class G, class C>
    inline auto structured_grid<G, C>::length() const noexcept
        -> length_type
    {
        return derived_grid().m_length;
    }

    /**
     * Returns a constant reference to the array of status at grid nodes.
     */
    template <class G, class C>
    inline auto structured_grid<G, C>::status_at_nodes() const
        -> const node_status_type&
    {
        return derived_grid().m_status_at_nodes;
    }

    /**
     * Returns a constant reference to the neighbors count of a given grid node.
     * 
     * @param idx Index of the grid node.
     */
    template <class G, class C>
    inline auto structured_grid<G, C>::neighbors_count(const size_type& idx) const
        -> const neighbors_count_type&
    {
        return derived_grid().neighbors_count(idx);
    }

    template <class G, class C>
    inline auto structured_grid<G, C>::neighbors_indices(const size_type& idx)
        -> neighbors_indices_type
    {
        neighbors_indices_type indices = xt::adapt(neighbors_indices_impl(idx));
        auto view = xt::view(indices, xt::range(0, neighbors_count(idx)));

        return view;
    }

    template <class G, class C>
    inline auto structured_grid<G, C>::neighbors_distances(const size_type& idx) const
        -> neighbors_distances_type
    {
        neighbors_distances_type distances = xt::adapt(neighbors_distances_impl(idx));
        auto view = xt::view(distances, xt::range(0, neighbors_count(idx)));

        return view;
    }

    /**
     * Iterate over the neighbors indices of a given grid node.
     *
     * Returns a fixed size std::array for performance considerations,
     * always use this method with the `neighbors_count` one.
     * 
     * Follows looped boundary conditions, if any.
     *
     * @param idx Index of the grid node.
     * @return Reference to the array of the neighbors indices of that grid node.
     */
    template <class G, class C>
    inline auto structured_grid<G, C>::neighbors_indices_impl(const size_type& idx)
        -> const neighbors_indices_impl_type&
    {
        if (m_neighbors_indices_cache.has(idx))
        {
            neighbors_indices_impl_type& n_indices = m_neighbors_indices_cache.get(idx);
            return n_indices;
        } else
        {
            neighbors_indices_impl_type& n_indices = m_neighbors_indices_cache.get_storage(idx);
            derived_grid().neighbors_indices_impl(n_indices, idx);
            return n_indices;
        } 
    }

    template <class G, class C>
    inline auto structured_grid<G, C>::neighbors_distances_impl(const size_type& idx) const
        -> const neighbors_distances_impl_type&
    {
        return derived_grid().neighbors_distances_impl(idx);
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
     * @return The vector of the neighbors of that grid node.
     */
    template <class G, class C>
    inline auto structured_grid<G, C>::neighbors(const size_type& idx)
        -> neighbors_type
    {
        neighbors_type node_neighbors;
        neighbors(idx, node_neighbors);

        return node_neighbors;
    }

    /**
     * Iterate over the neighbors indices of a given grid node.
     *
     * Follows looped boundary conditions, if any.
     *
     * @param idx Index of the grid node.
     * @param neighbors_indices Reference to the vector to be filled with the neighbors 
     *                          indices of that grid node.
     */
    template <class G, class C>
    inline void structured_grid<G, C>::neighbors_indices(const size_type& idx, neighbors_indices_type& neighbors_indices)
    {    
        const auto& n_count = neighbors_count(idx);
        const auto& n_indices = neighbors_indices_impl(idx);

        if (neighbors_indices.size() != n_count)
        {
            neighbors_indices.resize({n_count});
        }
        
        for (neighbors_count_type i=0; i<n_count; ++i)
        {   
            neighbors_indices[i] = n_indices[i];
        }
    }

    /**
     * Iterate over the neighbors of a given grid node.
     *
     * Follows looped boundary conditions, if any.
     *
     * @param idx Index of the grid node.
     * @param neighbors Reference to the vector to be filled with the neighbors of that grid node.
     */
    template <class G, class C>
    inline void structured_grid<G, C>::neighbors(const size_type& idx, neighbors_type& neighbors)
    {    
        size_type n_idx;
        const auto& n_count = neighbors_count(idx);
        const auto& n_indices = neighbors_indices_impl(idx);
        const auto& n_distances = neighbors_distances_impl(idx);

        if (neighbors.size() != n_count)
        {
            neighbors.resize({n_count});
        }
        
        for (neighbors_count_type i=0; i<n_count; ++i)
        {   
            n_idx = n_indices[i];
            neighbors[i] = neighbor({n_idx, n_distances[i], status_at_nodes()[n_idx]});
        }
    }
    //@}
}

#endif