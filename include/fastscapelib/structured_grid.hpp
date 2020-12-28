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
         * Provides a cache for neighbors queries of a particular node of the grid.
         */
        template <class G, unsigned int N>
        class neighbors_cache
        {
            using neighbors_type = xt::xtensor<neighbor, 1>;
            using neighbors_shape_type = typename neighbors_type::shape_type;
            using neighbors_indices_type = std::array<std::size_t, N>;
            using all_neighbors_ind_type = xt::xtensor<neighbors_indices_type, 1>;
            using all_neighbors_shape_type = typename all_neighbors_ind_type::shape_type;

            const G& m_grid;
            
            all_neighbors_ind_type m_all_neighbors_indices;

            neighbors_type m_neighbors = neighbors_type(neighbors_shape_type({N}));

        public:

            neighbors_cache(const G& grid, std::size_t size)
                : m_grid(grid), m_all_neighbors_indices(all_neighbors_shape_type({size}))
            {
                for (std::size_t i=0; i<size; ++i)
                {
                    m_all_neighbors_indices[i].fill(std::numeric_limits<std::size_t>::max());
                }
            }

            const neighbors_type& get_neighbors(const std::size_t idx)
            {
                using neighbors_distances_type = typename G::neighbors_distances_type;
                using neighbors_count_type = typename G::neighbors_count_type;
            
                neighbors_indices_type& neighbors_indices = m_all_neighbors_indices[idx];
                if(neighbors_indices[0] == std::numeric_limits<std::size_t>::max())
                {
                    m_grid.neighbors_indices_impl(neighbors_indices, idx);
                };
                
                std::size_t n_idx;
                neighbors_count_type neighbors_count = m_grid.neighbors_count(idx);
                const neighbors_distances_type& neighbors_distances = m_grid.neighbors_distance_impl(idx);

                if (neighbors_count != m_neighbors.size())
                {
                    m_neighbors.resize({neighbors_count});
                }
                
                for (neighbors_count_type i=0; i<neighbors_count; ++i)
                {   
                    n_idx = neighbors_indices[i];
                    m_neighbors[i] = neighbor({n_idx, neighbors_distances[i], m_grid.m_status_at_nodes[n_idx]});
                }

                return m_neighbors;
            }

            void warm_up() 
            {
                for (std::size_t i=0; i<m_grid.size(); ++i)
                {
                    neighbors_indices_type& all_indices = m_all_neighbors_indices[i];
                    m_grid.neighbors_indices_impl(all_indices, i);
                }
            };

            void reset()
            {
                for (std::size_t i=0; i<m_grid.size(); ++i)
                {
                    m_all_neighbors_indices[i].fill(std::numeric_limits<std::size_t>::max());
                }
            }

            std::size_t size() 
            {
                return m_all_neighbors_indices.size();
            };

            std::size_t used() 
            {
                std::size_t count = 0;
                
                for (std::size_t i=0; i < m_grid.size(); ++i)
                {
                    if (m_all_neighbors_indices[i][0] != std::numeric_limits<std::size_t>::max()) { count += 1; }
                }

                return count; 
            };
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
     * @tparam N Maximum neighbors count of a node.
     * @tparam C Neighbors cache type.
     */
    template <class G, unsigned int N, template<class, unsigned int> class C = detail::neighbors_cache>
    class structured_grid
    {
    public:

        using derived_grid_type = G;
        using inner_types = grid_inner_types<derived_grid_type>;

        using size_type = typename inner_types::size_type;
        using shape_type = typename inner_types::shape_type;
        using spacing_type = typename inner_types::spacing_type;
        
        using neighbors_indices_type = typename std::array<std::size_t, N>;
        using neighbors_distances_type = typename std::array<double, N>;
        using neighbors_count_type = typename inner_types::neighbors_count_type;
        using neighbors_type = typename xt::xtensor<neighbor, 1>;
        
        // using neighbors_distance_type = typename std::array<std::size_t, 8>;
        using neighbors_distance_type = xt::xtensor<double, 1>;

        //using neighbors_type = std::array<std::size_t, 8>;

        using spacing_t = std::conditional_t<std::is_arithmetic<spacing_type>::value,
                                            spacing_type,
                                            const spacing_type&>;

        using node_status_type = typename inner_types::node_status_type;
        using neighbors_cache_type = C<G, N>;
        
        //shape_type shape() const noexcept;
        size_type size() const noexcept;
        spacing_t spacing() const noexcept;
        const node_status_type& status_at_nodes() const;

        neighbors_count_type neighbors_count(const size_type& idx) const;

        const neighbors_type& neighbors(std::size_t idx);

        neighbors_cache_type neighbors_cache() { return m_neighbors_cache; };

    protected:

        structured_grid(std::size_t size)
            : m_neighbors_cache(neighbors_cache_type(derived_grid(), size))
        {};

        ~structured_grid() = default;

        neighbors_cache_type m_neighbors_cache;

        const derived_grid_type& derived_grid() const & noexcept;

        // friend C;
    };

    template <class G, unsigned int N, template<class, unsigned int> class C>
    inline auto structured_grid<G, N, C>::derived_grid() const & noexcept -> const derived_grid_type&
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
    template <class G, unsigned int N, template<class, unsigned int> class C >
    inline auto structured_grid<G, N, C>::size() const noexcept -> size_type
    {
        return derived_grid().m_size;
    }

    /**
     * Returns the (uniform) spacing between two adjacent grid nodes.
     *
     * Returns a copy of the value for 1-d grids or a constant reference otherwise.
     */
    template <class G, unsigned int N, template<class, unsigned int> class C >
    inline auto structured_grid<G, N, C>::spacing() const noexcept -> spacing_t
    {
        return derived_grid().m_spacing;
    }

    /**
     * Returns a constant reference to the array of status at grid nodes.
     */
    template <class G, unsigned int N, template<class, unsigned int> class C >
    inline auto structured_grid<G, N, C>::status_at_nodes() const -> const node_status_type&
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
    template <class G, unsigned int N, template<class, unsigned int> class C >
    inline auto structured_grid<G, N, C>::neighbors(std::size_t idx) -> const neighbors_type&
    {
        return m_neighbors_cache.get_neighbors(idx);
    }
    //@}
}