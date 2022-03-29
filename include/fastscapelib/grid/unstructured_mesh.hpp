#ifndef FASTSCAPELIB_GRID_UNSTRUCTURED_MESH_H
#define FASTSCAPELIB_GRID_UNSTRUCTURED_MESH_H

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{

    template <class S, class C>
    class unstructured_mesh_xt;

    template <class S, class C>
    struct grid_inner_types<unstructured_mesh_xt<S, C>>
    {
        static constexpr bool is_structured = false;
        static constexpr bool is_uniform = false;
        static constexpr uint8_t max_neighbors = 50;
        static constexpr std::size_t xt_ndims = 1;

        using xt_selector = S;
        using xt_type = xt_tensor_t<xt_selector, int, xt_ndims>;
        using field_data_type = double;

        using size_type = typename xt_type::size_type;
        using shape_type = typename xt_type::shape_type;
        using distance_type = double;

        using neighbors_cache_type = C;
        using neighbors_count_type = std::uint8_t;
        using neighbors_distances_impl_type = typename std::array<distance_type, max_neighbors>;

        using node_status_type = xt_tensor_t<xt_selector, node_status, xt_ndims>;
    };

    /**
     * @class unstructured_mesh_xt
     * @brief 2-dimensional unstructured mesh.
     *
     * @tparam S xtensor container selector for data array members.
     * @tparam C Grid neighbors cache type.
     */
    template <class S, class C = neighbors_no_cache<50>>
    class unstructured_mesh_xt : public grid<unstructured_mesh_xt<S, C>>
    {
    public:
        using self_type = unstructured_mesh_xt<S, C>;
        using base_type = grid<self_type>;
        using inner_types = grid_inner_types<self_type>;

        using xt_selector = typename inner_types::xt_selector;
        static constexpr std::size_t xt_ndims = inner_types::xt_ndims;
        static constexpr std::uint8_t max_neighbors()
        {
            return inner_types::max_neighbors;
        };

        using size_type = typename inner_types::size_type;
        using shape_type = typename inner_types::shape_type;
        using distance_type = typename inner_types::distance_type;

        using neighbors_type = typename base_type::neighbors_type;
        using neighbors_count_type = typename inner_types::neighbors_count_type;
        using neighbors_indices_type = typename base_type::neighbors_indices_type;
        using neighbors_distances_type = typename base_type::neighbors_distances_type;

        using node_status_type = typename inner_types::node_status_type;

        unstructured_mesh_xt()
            : base_type(0){};

    protected:
        using neighbors_distances_impl_type = typename inner_types::neighbors_distances_impl_type;
        using neighbors_indices_impl_type = typename base_type::neighbors_indices_impl_type;

        shape_type m_shape;
        size_type m_size;
        distance_type m_node_area;

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
