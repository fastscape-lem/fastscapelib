#ifndef FASTSCAPELIB_GRID_UNSTRUCTURED_MESH_H
#define FASTSCAPELIB_GRID_UNSTRUCTURED_MESH_H

#include "xtensor/xindex_view.hpp"

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"

#include <cmath>


namespace fastscapelib
{

    template <class S, unsigned int N>
    class unstructured_mesh_xt;

    template <class S, unsigned int N>
    struct grid_inner_types<unstructured_mesh_xt<S, N>>
    {
        static constexpr bool is_structured = false;
        static constexpr bool is_uniform = false;

        using grid_data_type = double;

        using xt_selector = S;
        static constexpr std::size_t xt_ndims = 1;

        static constexpr uint8_t n_neighbors_max = N;
        using neighbors_cache_type = neighbors_no_cache<0>;
        using neighbors_count_type = std::uint8_t;
    };

    /**
     * @class unstructured_mesh_xt
     * @brief 2-dimensional unstructured mesh.
     *
     * @tparam S xtensor container selector for data array members.
     * @tparam N Max number of grid node neighbors.
     */
    template <class S, unsigned int N = 30>
    class unstructured_mesh_xt : public grid<unstructured_mesh_xt<S, N>>
    {
    public:
        using self_type = unstructured_mesh_xt<S, N>;
        using base_type = grid<self_type>;

        using grid_data_type = typename base_type::grid_data_type;

        using xt_selector = typename base_type::xt_selector;
        using size_type = typename base_type::size_type;
        using shape_type = typename base_type::shape_type;

        using points_type = xt_tensor_t<xt_selector, grid_data_type, 2>;
        using indices_type = xt_tensor_t<xt_selector, size_type, 1>;
        using areas_type = xt_tensor_t<xt_selector, grid_data_type, 1>;

        using neighbors_type = typename base_type::neighbors_type;
        using neighbors_count_type = typename base_type::neighbors_count_type;
        using neighbors_indices_type = typename base_type::neighbors_indices_type;
        using neighbors_distances_type = typename base_type::neighbors_distances_type;

        using node_status_type = typename base_type::node_status_type;

        void set_status_at_nodes(const std::vector<node>& status_at_nodes);

        const neighbors_count_type neighbors_count(const size_type& idx) const;

        unstructured_mesh_xt(const points_type& points,
                             const indices_type& neighbors_indices_ptr,
                             const indices_type& neighbors_indices,
                             const indices_type& convex_hull_indices,
                             const areas_type& areas,
                             const std::vector<node>& status_at_nodes = {});

    protected:
        using neighbors_distances_impl_type = typename base_type::neighbors_distances_impl_type;
        using neighbors_indices_impl_type = typename base_type::neighbors_indices_impl_type;

        shape_type m_shape;
        size_type m_size;
        neighbors_distances_type m_distances; 
        grid_data_type m_node_area;

        points_type m_points;
        indices_type m_neighbors_indices_ptr;
        indices_type m_neighbors_indices;
        indices_type m_convex_hull_indices;
        areas_type m_areas;

        node_status_type m_status_at_nodes;

        void neighbors_indices_impl(neighbors_indices_impl_type& neighbors,
                                    const size_type& idx) const;

        inline const neighbors_distances_impl_type& neighbors_distances_impl(
            const size_type& idx) const;

        friend class grid<self_type>;

        // void neighbors_distances_impl(neighbors_distances_impl_type& neighbors_distances, 
        //                             const size_type& idx) const;

        neighbors_distances_type& neighbors_distances[m_size];

        for (const auto x : points):
            idx = neighbors[x] //uses the neighbors derived from indices (i.e. from scipy neighbors = indices[indptr[x]:indptr[x+1]])
            for (const auto y : neighbors);
                neighbors_distances = sqrt((pow((pointstest[neighbors[y]][0])-(points[x][0])))+(pow((pointstest[neighbors[y]][1])-(points[x][1]))));


    };


    template <class S, unsigned int N>
    unstructured_mesh_xt<S, N>::unstructured_mesh_xt(const points_type& points,
                                                     const indices_type& neighbors_indices_ptr,
                                                     const indices_type& neighbors_indices,
                                                     const indices_type& convex_hull_indices,
                                                     const areas_type& areas,
                                                     const std::vector<node>& status_at_nodes)
        // no neighbors cache -> base type argument value doesn't matter
        : base_type(0)
        , m_points(points)
        , m_neighbors_indices_ptr(neighbors_indices_ptr)
        , m_neighbors_indices(neighbors_indices)
        , m_convex_hull_indices(convex_hull_indices)
        , m_areas(areas)
    {
        // TODO: check that all array shapes are consistent?

        m_size = points.shape()[0];
        m_shape = { static_cast<typename shape_type::value_type>(m_size) };
        set_status_at_nodes(status_at_nodes);
        // TODO: pre-compute neighbor distances
        m_neighbors_distances = neighbors_distances;
    }

    template <class S, unsigned int N>
    void unstructured_mesh_xt<S, N>::set_status_at_nodes(const std::vector<node>& status_at_nodes)
    {
        node_status_type temp_status_at_nodes(m_shape, node_status::core);

        if (status_at_nodes.size() > 0)
        {
            for (const node& n : status_at_nodes)
            {
                if (n.status == node_status::looped_boundary)
                {
                    throw std::invalid_argument("node_status::looped_boundary is not allowed in "
                                                "unstructured meshes");
                }

                temp_status_at_nodes.at(n.idx) = n.status;
            }
        }
        else
        {
            // if no status at node is given, set fixed value boundaries for all nodes
            // forming the convex hull
            auto status_at_qhull_view = xt::index_view(temp_status_at_nodes, m_convex_hull_indices);
            status_at_qhull_view = node_status::fixed_value_boundary;
        }

        m_status_at_nodes = temp_status_at_nodes;
    }

    template <class S, unsigned int N>
    auto unstructured_mesh_xt<S, N>::neighbors_count(const size_type& idx) const
        -> const neighbors_count_type
    {
        size_type start_idx = m_neighbors_indices_ptr[idx];
        size_type stop_idx = m_neighbors_indices_ptr[idx + 1];

        return static_cast<neighbors_count_type>(stop_idx - start_idx);
    }


    template <class S, unsigned int N>
    void unstructured_mesh_xt<S, N>::neighbors_indices_impl(neighbors_indices_impl_type& neighbors,
                                                            const size_type& idx) const
    {
        int start_idx = m_neighbors_indices_ptr[idx];
        int stop_idx = m_neighbors_indices_ptr[idx + 1];

        for (auto i = 0; i == (stop_idx - start_idx);)
        {
            neighbors[i] = m_neighbors_indices[start_idx + i];
        }
    }

    template <class S, unsigned int N>
    auto unstructured_mesh_xt<S, N>::neighbors_distances_impl(const size_type& idx) const
        -> const neighbors_distances_impl_type&
    {
        return neighbors_distances_impl_type();  // for this look at indices, this should hold an
                                                 // array of pre-computed distances
    }


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
