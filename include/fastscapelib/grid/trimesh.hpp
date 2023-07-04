#ifndef FASTSCAPELIB_GRID_TRIMESH_H
#define FASTSCAPELIB_GRID_TRIMESH_H

#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "xtensor/xindex_view.hpp"

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{

    namespace detail
    {

        /**
         * Used to extract unique edges in the mesh (or count their occurence).
         *
         * This hash function yields the same output for simple permutations,
         * which is what we want since the pairs of node indices (2, 5) and
         * (5, 2) both refer to the same edge.
         */
        template <class T>
        struct tri_edge_hash
        {
            std::size_t operator()(const std::pair<T, T>& p) const
            {
                auto h1 = std::hash<T>()(p.first);
                auto h2 = std::hash<T>()(p.second);

                return h1 ^ h2;
            }
        };

        /**
         * Used to extract unique edges in the mesh (or count their occurence).
         *
         * This comparator returns true for simple permuations. It is needed in
         * addition to ``edge_hash`` in order to build a map of unique edges
         * ignoring permutations.
         */
        template <class T>
        struct tri_edge_equal
        {
            using pair_type = std::pair<T, T>;

            bool operator()(const pair_type& p1, const pair_type& p2) const
            {
                if (p1.first == p2.second && p1.second == p2.first)
                {
                    return true;
                }
                else
                {
                    return p1.first == p2.first && p1.second == p2.second;
                }
            }
        };

        template <class T>
        using tri_edge_type = std::pair<T, T>;

        template <class T>
        using tri_edge_map = std::
            unordered_map<tri_edge_type<T>, std::size_t, tri_edge_hash<T>, tri_edge_equal<T>>;
    }


    template <class S, unsigned int N>
    class trimesh_xt;

    /**
     * Unstructured mesh specialized types.
     */
    template <class S, unsigned int N>
    struct grid_inner_types<trimesh_xt<S, N>>
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
     * @brief 2-dimensional unstructured mesh.
     *
     * Fastscapelib grid adapter for a 2-d triangulated mesh. This class
     * requires an input mesh (it doesn't provide any meshing capability).
     *
     * @tparam S The xtensor container selector for data array members.
     * @tparam N The maximum number of node neighbors.
     */
    template <class S, unsigned int N = 20>
    class trimesh_xt : public grid<trimesh_xt<S, N>>
    {
    public:
        using self_type = trimesh_xt<S, N>;
        using base_type = grid<self_type>;

        using grid_data_type = typename base_type::grid_data_type;

        using xt_selector = typename base_type::xt_selector;
        using size_type = typename base_type::size_type;
        using shape_type = typename base_type::shape_type;

        using points_type = xt_tensor_t<xt_selector, grid_data_type, 2>;
        using triangles_type = xt_tensor_t<xt_selector, size_type, 2>;
        using indices_type = xt_tensor_t<xt_selector, size_type, 1>;
        using areas_type = xt_tensor_t<xt_selector, grid_data_type, 1>;

        using neighbors_type = typename base_type::neighbors_type;
        using neighbors_count_type = typename base_type::neighbors_count_type;
        using neighbors_indices_type = typename base_type::neighbors_indices_type;
        using neighbors_distances_type = typename base_type::neighbors_distances_type;

        using node_status_type = typename base_type::node_status_type;

        const neighbors_count_type& neighbors_count(const size_type& idx) const;

        trimesh_xt(const points_type& points,
                   const triangles_type& triangles,
                   const areas_type& areas,
                   const std::vector<node>& nodes_status = {});

    protected:
        using neighbors_distances_impl_type = typename base_type::neighbors_distances_impl_type;
        using neighbors_indices_impl_type = typename base_type::neighbors_indices_impl_type;

        shape_type m_shape;
        size_type m_size;

        points_type m_nodes_points;
        areas_type m_nodes_areas;
        std::unordered_set<size_type> m_boundary_nodes;
        node_status_type m_nodes_status;

        std::vector<neighbors_indices_impl_type> m_neighbors_indices;
        std::vector<neighbors_distances_impl_type> m_neighbors_distances;
        std::vector<neighbors_count_type> m_neighbors_counts;

        void set_nodes_status(const std::vector<node>& nodes_status);

        inline areas_type nodes_areas_impl() const;
        inline grid_data_type nodes_areas_impl(const size_type& idx) const noexcept;

        void neighbors_indices_impl(neighbors_indices_impl_type& neighbors,
                                    const size_type& idx) const;

        inline const neighbors_distances_impl_type& neighbors_distances_impl(
            const size_type& idx) const;

        friend class grid<self_type>;
    };


    /**
     * @name Constructors
     */
    //@{
    /**
     * Creates a new mesh.
     *
     * @param points The mesh node x,y coordinates (array of shape [N, 2]).
     * @param triangles The node indices of the triangles (array of shape [K, 3]).
     * @param areas The area of the cells centered to each mesh node (array of shape [N]).
     * @param nodes_status Manually define the status at any node on the mesh.
     *
     * If ``nodes_status`` is empty, a "fixed value" status is set for all
     * boundary nodes (i.e., the end-points of all the edges that are not shared
     * by more than one triangle).
     */
    template <class S, unsigned int N>
    trimesh_xt<S, N>::trimesh_xt(const points_type& points,
                                 const triangles_type& triangles,
                                 const areas_type& areas,
                                 const std::vector<node>& nodes_status)
        : base_type(0)
        , m_nodes_points(points)
        , m_nodes_areas(areas)
    {
        // TODO: sanity checks, e.g., all array shapes are consistent

        m_size = points.shape()[0];
        m_shape = { static_cast<typename shape_type::value_type>(m_size) };

        // extract and count triangle edges

        using edge_type = detail::tri_edge_type<size_type>;
        using edge_map = detail::tri_edge_map<size_type>;

        edge_map edges_count;
        const std::array<std::array<size_type, 2>, 3> tri_local_indices{
            { { 1, 2 }, { 2, 0 }, { 0, 1 } }
        };

        size_type n_triangles = triangles.shape()[0];

        for (size_type i = 0; i < n_triangles; i++)
        {
            for (const auto& edge_idx : tri_local_indices)
            {
                const edge_type key(triangles(i, edge_idx[0]), triangles(i, edge_idx[1]));

                auto result = edges_count.insert({ key, 1 });
                // increment edge count if already inserted
                if (!result.second)
                {
                    result.first->second += 1;
                }
            }
        }

        // fill node neighbor data and find boundary nodes

        m_boundary_nodes.clear();
        m_neighbors_counts.assign(m_size, 0);
        m_neighbors_indices.resize(m_size);
        m_neighbors_distances.resize(m_size);

        for (const auto& edge : edges_count)
        {
            const edge_type& edge_points = edge.first;
            size_type count = edge.second;

            if (count == 1)
            {
                m_boundary_nodes.insert(edge_points.first);
                m_boundary_nodes.insert(edge_points.second);
            }

            m_neighbors_counts[edge_points.first]++;
            m_neighbors_counts[edge_points.second]++;
            m_neighbors_indices[edge_points.first].push_back(edge_points.second);
            m_neighbors_indices[edge_points.second].push_back(edge_points.first);

            const auto x1 = m_nodes_points(edge_points.first, 0);
            const auto y1 = m_nodes_points(edge_points.first, 1);
            const auto x2 = m_nodes_points(edge_points.second, 0);
            const auto y2 = m_nodes_points(edge_points.second, 1);
            auto distance = std::sqrt(((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)));
            m_neighbors_distances[edge_points.first].push_back(distance);
            m_neighbors_distances[edge_points.second].push_back(distance);
        }

        set_nodes_status(nodes_status);
    }
    //@}

    template <class S, unsigned int N>
    void trimesh_xt<S, N>::set_nodes_status(const std::vector<node>& nodes_status)
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
            // if no status at node is given, set fixed value boundaries for all boundary nodes
            for (const size_type& idx : m_boundary_nodes)
            {
                temp_nodes_status[idx] = node_status::fixed_value;
            }
        }

        m_nodes_status = temp_nodes_status;
    }

    /**
     * @name Neighbor methods
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
    auto trimesh_xt<S, N>::neighbors_count(const size_type& idx) const
        -> const neighbors_count_type&
    {
        return m_neighbors_counts[idx];
    }
    //@}

    template <class S, unsigned int N>
    inline auto trimesh_xt<S, N>::nodes_areas_impl() const -> areas_type
    {
        return m_nodes_areas;
    }

    template <class S, unsigned int N>
    inline auto trimesh_xt<S, N>::nodes_areas_impl(const size_type& idx) const noexcept
        -> grid_data_type
    {
        return m_nodes_areas(idx);
    }

    template <class S, unsigned int N>
    void trimesh_xt<S, N>::neighbors_indices_impl(neighbors_indices_impl_type& neighbors,
                                                  const size_type& idx) const
    {
        const auto& size = m_neighbors_counts[idx];
        neighbors.resize(size);

        for (size_type i = 0; i < size; i++)
        {
            neighbors[i] = m_neighbors_indices[idx][i];
        }
    }

    template <class S, unsigned int N>
    auto trimesh_xt<S, N>::neighbors_distances_impl(const size_type& idx) const
        -> const neighbors_distances_impl_type&
    {
        return m_neighbors_distances[idx];
    }


    /**
     * @typedef trimesh
     *
     * \rst
     * Alias template on ``trimesh_xt`` with :cpp:type:`xt::xtensor`
     * used as array container type for data members.
     *
     * This is mainly for convenience when using in C++ applications.
     * \endrst
     */
    using trimesh = trimesh_xt<xt_selector>;
}

#endif  // FASTSCAPELIB_GRID_TRIMESH_H
