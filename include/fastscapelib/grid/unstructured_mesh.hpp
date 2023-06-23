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

    /**
     * Unstructured mesh specialized types.
     */
    template <class S, unsigned int N>
    struct grid_inner_types<unstructured_mesh_xt<S, N>>
    {
        static constexpr bool is_structured = false;
        static constexpr bool is_uniform = false;

        using grid_data_type = double;

        using xt_selector = S;
        static constexpr std::size_t xt_ndims = 1;

        static constexpr uint8_t n_neighbors_max = N;
        using neighbors_cache_type = neighbors_no_cache<N>;
        using neighbors_count_type = std::uint8_t;
    };

    /**
     * @brief 2-dimensional unstructured mesh.
     *
     * Fastscapelib grid adapter for a 2-d triangulated mesh. This class
     * requires an input mesh (it doesn't provide any meshing capability).
     *
     * @tparam S The xtensor container selector for data array members.
     * @tparam N The maximum number of grid node neighbors.
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

        void set_nodes_status(const std::vector<node>& nodes_status);

        inline grid_data_type nodes_areas(const size_type& idx) const noexcept;
        const neighbors_count_type& neighbors_count(const size_type& idx) const;

        unstructured_mesh_xt(const points_type& points,
                             const indices_type& neighbors_indices_ptr,
                             const indices_type& neighbors_indices,
                             const indices_type& convex_hull_indices,
                             const areas_type& areas,
                             const std::vector<node>& nodes_status = {});

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

        node_status_type m_nodes_status;

        std::vector<neighbors_distances_impl_type> m_neighbors_distances;
        std::vector<neighbors_count_type> m_neighbors_counts;

        void neighbors_indices_impl(neighbors_indices_impl_type& neighbors,
                                    const size_type& idx) const;

        inline const neighbors_distances_impl_type& neighbors_distances_impl(
            const size_type& idx) const;

        void compute_neighbors_counts();
        void compute_neighbors_distances();

        friend class grid<self_type>;
    };


    /**
     * @name Constructors
     */
    //@{
    /**
     * Creates a new mesh.
     *
     * @param points The mesh node x,y coordinates (expects an array of shape [N, 2]).
     * @param neighbors_indices_ptr The lookup offsets of the neighbors of each mesh node
     *                              (expects an array of shape [N+1]).
     * @param neighbors_indices The node neighbor indices (flattened array)
     * @param convex_hull_indices The indices of the boundary nodes.
     * @param areas The area of the cells centered to each mesh node (array of shape [N]).
     * @param nodes_status Manually define the status at any node on the mesh.
     *
     * If ``nodes_status`` is empty, a "fixed value" status is set for all
     * boundary nodes.
     */
    template <class S, unsigned int N>
    unstructured_mesh_xt<S, N>::unstructured_mesh_xt(const points_type& points,
                                                     const indices_type& neighbors_indices_ptr,
                                                     const indices_type& neighbors_indices,
                                                     const indices_type& convex_hull_indices,
                                                     const areas_type& areas,
                                                     const std::vector<node>& nodes_status)
        // no neighbors cache -> base type argument value doesn't matter
        : base_type(0)
        , m_points(points)
        , m_neighbors_indices_ptr(neighbors_indices_ptr)
        , m_neighbors_indices(neighbors_indices)
        , m_convex_hull_indices(convex_hull_indices)
        , m_areas(areas)
    {
        // TODO: sanity checks, e.g., all array shapes are consistent

        m_size = points.shape()[0];
        m_shape = { static_cast<typename shape_type::value_type>(m_size) };
        set_nodes_status(nodes_status);

        // counts must be called before distances, both after setting m_size
        compute_neighbors_counts();
        compute_neighbors_distances();
    }
    //@}

    // pre-compute and store the number of neighbors wastes memory for such trivial operation
    // but it is needed since the `neighbors_count` public method returns a const reference
    template <class S, unsigned int N>
    void unstructured_mesh_xt<S, N>::compute_neighbors_counts()
    {
        m_neighbors_counts.clear();
        m_neighbors_counts.resize(m_size);

        for (size_t i = 0; i < m_size; i++)
        {
            size_type start_idx = m_neighbors_indices_ptr[i];
            size_type stop_idx = m_neighbors_indices_ptr[i + 1];

            m_neighbors_counts[i] = static_cast<neighbors_count_type>(stop_idx - start_idx);
        }
    }

    template <class S, unsigned int N>
    void unstructured_mesh_xt<S, N>::compute_neighbors_distances()
    {
        m_neighbors_distances.clear();
        m_neighbors_distances.resize(m_size);

        for (size_t i = 0; i < m_size; i++)
        {
            const auto ix = m_points(i, 0);
            const auto iy = m_points(i, 1);

            size_type start_idx = m_neighbors_indices_ptr[i];

            neighbors_distances_impl_type nb_distances;

            for (size_type inb = 0; inb < m_neighbors_counts[i]; inb++)
            {
                const auto nb_idx = m_neighbors_indices[start_idx + inb];
                const auto nbx = m_points(nb_idx, 0);
                const auto nby = m_points(nb_idx, 1);
                nb_distances[inb] = std::sqrt(((nbx - ix) * (nbx - ix) + (nby - iy) * (nby - iy)));
            }

            m_neighbors_distances[i] = nb_distances;
        }
    }

    template <class S, unsigned int N>
    void unstructured_mesh_xt<S, N>::set_nodes_status(const std::vector<node>& nodes_status)
    {
        node_status_type temp_nodes_status(m_shape, node_status::core);

        if (nodes_status.size() > 0)
        {
            for (const node& n : nodes_status)
            {
                if (n.status == node_status::looped)
                {
                    throw std::invalid_argument("node_status::looped is not allowed in "
                                                "unstructured meshes");
                }

                temp_nodes_status.at(n.idx) = n.status;
            }
        }
        else
        {
            // if no status at node is given, set fixed value boundaries for all nodes
            // forming the convex hull
            auto status_at_qhull_view = xt::index_view(temp_nodes_status, m_convex_hull_indices);
            status_at_qhull_view = node_status::fixed_value;
        }

        m_nodes_status = temp_nodes_status;
    }

    /**
     * @name Grid properties
     */
    //@{
    template <class S, unsigned int N>
    inline auto unstructured_mesh_xt<S, N>::nodes_areas(const size_type& idx) const noexcept
        -> grid_data_type
    {
        return m_areas[idx];
    }
    //@}

    /**
     * @name Grid topology
     */
    /**
     * Returns the number of neighbors of a given grid node.
     *
     * @param idx The grid node flat index.
     *
     * @see fastscapelib::grid<G>::neighbors_indices,
     *      fastscapelib::grid<G>::neighbors_distances,
     *      fastscapelib::grid<G>::neighbors
     */
    template <class S, unsigned int N>
    auto unstructured_mesh_xt<S, N>::neighbors_count(const size_type& idx) const
        -> const neighbors_count_type&
    {
        return m_neighbors_counts[idx];
    }
    //@}

    template <class S, unsigned int N>
    void unstructured_mesh_xt<S, N>::neighbors_indices_impl(neighbors_indices_impl_type& neighbors,
                                                            const size_type& idx) const
    {
        size_type start_idx = m_neighbors_indices_ptr[idx];
        size_type stop_idx = m_neighbors_indices_ptr[idx + 1];

        for (size_type i = 0; i < (stop_idx - start_idx); i++)
        {
            neighbors[i] = m_neighbors_indices[start_idx + i];
        }
    }

    template <class S, unsigned int N>
    auto unstructured_mesh_xt<S, N>::neighbors_distances_impl(const size_type& idx) const
        -> const neighbors_distances_impl_type&
    {
        return m_neighbors_distances[idx];
    }


    /**
     * @typedef unstructured_mesh
     *
     * \rst
     * Alias template on ``unstructured_mesh_xt`` with :cpp:type:`xt::xtensor`
     * used as array container type for data members.
     *
     * This is mainly for convenience when using in C++ applications.
     * \endrst
     */
    using unstructured_mesh = unstructured_mesh_xt<xt_selector>;
}

#endif  // FASTSCAPELIB_GRID_UNSTRUCTURED_MESH_H
