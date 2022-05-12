#ifndef FASTSCAPELIB_GRID_UNSTRUCTURED_MESH_H
#define FASTSCAPELIB_GRID_UNSTRUCTURED_MESH_H

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{

    template <class S, unsigned int N, class C>
    class unstructured_mesh_xt;

    template <class S, unsigned int N, class C>
    struct grid_inner_types<unstructured_mesh_xt<S, N, C>>
    {
        static constexpr bool is_structured = false;
        static constexpr bool is_uniform = false;

        using grid_data_type = double;

        using xt_selector = S;
        static constexpr std::size_t xt_ndims = 1;

        static constexpr uint8_t max_neighbors = N;
        using neighbors_cache_type = C;
        using neighbors_count_type = std::uint8_t;
    };

    /**
     * @class unstructured_mesh_xt
     * @brief 2-dimensional unstructured mesh.
     *
     * @tparam S xtensor container selector for data array members.
     * @tparam N Max number of grid node neighbors.
     * @tparam C Grid neighbors cache type.
     */
    template <class S, unsigned int N = 50, class C = neighbors_no_cache<N>>
    class unstructured_mesh_xt : public grid<unstructured_mesh_xt<S, N, C>>
    {
    public:
        using self_type = unstructured_mesh_xt<S, N, C>;
        using base_type = grid<self_type>;

        using grid_data_type = typename base_type::grid_data_type;

        using xt_selector = typename base_type::xt_selector;
        using size_type = typename base_type::size_type;
        using shape_type = typename base_type::shape_type;

        using neighbors_type = typename base_type::neighbors_type;
        using neighbors_count_type = typename base_type::neighbors_count_type;
        using neighbors_indices_type = typename base_type::neighbors_indices_type;
        using neighbors_distances_type = typename base_type::neighbors_distances_type;

        using node_status_type = typename base_type::node_status_type;

        /**
         * Does it need a self type?
         * Why are the type classes defined in the protected area here but not in raster_grid.hpp?
         * is it ok that neighbors_indices is overloaded?
         */

        unstructured_mesh_xt(const shape_type& points,
                             const neighbors_count_type& neighbors_indices_ptr,
                             const neighbors_indices_type& neighbors_indices,
                             const neighbors_distances_type& convex_hull_indices,
                             const size_type& areas,
                             const node_status_type& status_at_nodes = {});


    protected:
        using neighbors_distances_impl_type = typename base_type::neighbors_distances_impl_type;
        using neighbors_indices_impl_type = typename base_type::neighbors_indices_impl_type;

        shape_type m_shape;
        size_type m_size;
        neighbors_distances_type m_spacing;
        grid_data_type m_node_area;

        std::array<double, grid_inner_types<self_type>::xt_ndims> m_points;
        std::array<double, grid_inner_types<self_type>::xt_ndims> m_neighbors_indices_ptr;
        std::array<double, grid_inner_types<self_type>::xt_ndims> m_neighbors_indices;
        std::array<double, grid_inner_types<self_type>::xt_ndims> m_convex_hull_indices;
        std::array<double, grid_inner_types<self_type>::xt_ndims> m_areas;

        // coded_ndistances_type m_neighbor_distances;

        node_status_type m_status_at_nodes;

        void set_status_at_nodes(const std::vector<node>& status_at_nodes);
        inline const neighbors_count_type& neighbors_count(const size_type& idx) const noexcept;

        void neighbors_indices_impl(neighbors_indices_impl_type& neighbors,
                                    const size_type& idx) const;

        const neighbors_distances_impl_type& neighbors_distances_impl(const size_type& idx) const;

        friend class grid<self_type>;
    };

    /**
     * @typedef unstructured_mesh
     * Alias template on unstructured_mesh_xt with ``xt::xtensor`` used
     * as array container type for data members.
     *
     * This is mainly for convenience when using in C++ applications.
     *
     */
    using unstructured_mesh = unstructured_mesh_xt<xt_selector>;
}

#endif  // FASTSCAPELIB_GRID_UNSTRUCTURED_MESH_H
