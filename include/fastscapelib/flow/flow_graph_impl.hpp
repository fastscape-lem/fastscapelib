#ifndef FASTSCAPELIB_FLOW_GRAPH_IMPL_H_
#define FASTSCAPELIB_FLOW_GRAPH_IMPL_H_

#include <array>
#include <iterator>
#include <limits>
#include <unordered_set>
#include <stack>
#include <vector>
#include <set>

#include "xtensor/xbroadcast.hpp"
#include "xtensor/xstrided_view.hpp"
#include "xtensor/xview.hpp"

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/utils/iterators.hpp"
#include "fastscapelib/utils/containers.hpp"


namespace fastscapelib
{

    namespace detail
    {

        /**
         * Flow operator implementation forward declaration.
         */
        template <class FG, class OP, class Tag>
        class flow_operator_impl;

        /**
         * Flow graph implementation.
         *
         * Implements the graph underlying data structures and API
         * (excepted for graph update that is implemented in flow operators).
         *
         * This template class has no definition, it must be specialized
         * using the Tag parameter.
         *
         * @tparam G The grid type
         * @tparam S The xtensor selector type
         * @tparam Tag The flow graph implementation Tag
         *
         */
        template <class G, class S, class Tag>
        class flow_graph_impl;
    }

    /**
     * Fixed array flow graph implementation.
     *
     * This implementation uses fixed shape arrays to store the graph
     * topology (node receivers, donors, etc.).
     *
     */
    struct flow_graph_fixed_array_tag
    {
    };


    namespace detail
    {

        template <class G, class S>
        class flow_graph_impl<G, S, flow_graph_fixed_array_tag>
        {
        public:
            using self_type = flow_graph_impl<G, S, flow_graph_fixed_array_tag>;
            using grid_type = G;
            using container_selector = S;
            using tag = flow_graph_fixed_array_tag;

            using size_type = typename grid_type::size_type;
            using data_type = typename grid_type::grid_data_type;
            using data_array_type = dynamic_shape_container_t<container_selector, data_type>;

            using donors_type = fixed_shape_container_t<container_selector, size_type, 2>;
            using donors_count_type = fixed_shape_container_t<container_selector, size_type, 1>;
            using receivers_type = donors_type;
            using receivers_count_type = donors_count_type;

            using receivers_distance_type
                = fixed_shape_container_t<container_selector, data_type, 2>;
            using receivers_weight_type = fixed_shape_container_t<container_selector, data_type, 2>;
            using nodes_indices_type = fixed_shape_container_t<container_selector, size_type, 1>;

            using basins_type = fixed_shape_container_t<container_selector, size_type, 1>;

            flow_graph_impl(grid_type& grid, bool single_flow = false)
                : m_single_flow(single_flow)
                , m_grid(grid)
            {
                size_type n_receivers_max = grid_type::n_neighbors_max();

                if (single_flow)
                {
                    n_receivers_max = 1;
                }

                using shape_type = std::array<size_type, 2>;
                const shape_type receivers_shape = { grid.size(), n_receivers_max };
                const shape_type donors_shape = { grid.size(), grid_type::n_neighbors_max() + 1 };

                m_receivers = xt::ones<size_type>(receivers_shape) * -1;
                m_receivers_count = xt::zeros<size_type>({ grid.size() });
                m_receivers_distance = xt::ones<data_type>(receivers_shape) * -1;
                m_receivers_weight = xt::zeros<data_type>(receivers_shape);

                m_donors = xt::ones<size_type>(donors_shape) * -1;
                m_donors_count = xt::zeros<size_type>({ grid.size() });

                // TODO: replace indices and levels for any traversal direction to save memory
                // footprint
                m_storage_indices = xt::arange<size_type>(0, grid.size(), 1);
                m_any_order_levels = nodes_indices_type({ 0, size() });

                m_dfs_indices = xt::ones<size_type>({ grid.size() }) * -1;
                m_bfs_indices = xt::ones<size_type>({ grid.size() }) * -1;
                m_bfs_levels = xt::ones<size_type>({ grid.size() + 1 }) * -1;

                // TODO: basins are not always needed (only init on-demand)
                m_basins = xt::empty<size_type>({ grid.size() });
            };

            bool single_flow() const
            {
                return m_single_flow;
            }

            G& grid() const
            {
                return m_grid;
            };

            size_type size() const
            {
                return m_grid.size();
            };

            const receivers_type& receivers() const
            {
                return m_receivers;
            };

            const receivers_count_type& receivers_count() const
            {
                return m_receivers_count;
            };

            const receivers_distance_type& receivers_distance() const
            {
                return m_receivers_distance;
            };

            const receivers_weight_type& receivers_weight() const
            {
                return m_receivers_weight;
            };

            const donors_type& donors() const
            {
                return m_donors;
            };

            const donors_count_type& donors_count() const
            {
                return m_donors_count;
            };

            void compute_donors();

            const nodes_indices_type& storage_indices() const
            {
                return m_storage_indices;
            };

            const nodes_indices_type& any_order_levels() const
            {
                return m_any_order_levels;
            };

            const nodes_indices_type& dfs_indices() const
            {
                return m_dfs_indices;
            };

            const nodes_indices_type& bfs_indices() const
            {
                return m_bfs_indices;
            };

            const nodes_indices_type& bfs_levels() const
            {
                return m_bfs_levels;
            };

            void compute_dfs_indices_bottomup();
            void compute_dfs_indices_topdown();
            void compute_bfs_indices_bottomup();

            /*
             * Should be used for graph traversal in an explicit direction.
             */
            inline stl_container_iterator_wrapper<nodes_indices_type> nodes_indices_bottomup() const
            {
                return m_dfs_indices;
            }

            const std::vector<size_type>& outlets() const
            {
                return m_outlets;
            }

            const std::vector<size_type>& pits();

            const basins_type& basins() const
            {
                return m_basins;
            };

            void compute_basins();

            template <class T>
            void accumulate(data_array_type& acc, T&& src) const;

            template <class T>
            data_array_type accumulate(T&& src) const;

            inline bool is_base_level(const size_type& idx) const
            {
                return bool(m_base_levels.count(idx));
            }

            const std::unordered_set<size_type>& base_levels() const
            {
                return m_base_levels;
            }

            template <class C>
            void set_base_levels(const C& levels)
            {
                m_base_levels.clear();
                m_base_levels.insert(levels.begin(), levels.end());
            }

            inline bool is_masked(const size_type& idx) const
            {
                return m_mask_initialized && m_mask.flat(idx);
            }

            const xt::xarray<bool>& mask() const
            {
                return m_mask;
            }

            template <class C>
            void set_mask(C&& mask)
            {
                m_mask = mask;
                m_mask_initialized = true;
            }

        private:
            bool m_single_flow;
            grid_type& m_grid;

            donors_type m_donors;
            donors_count_type m_donors_count;

            receivers_type m_receivers;
            receivers_count_type m_receivers_count;
            receivers_distance_type m_receivers_distance;
            receivers_weight_type m_receivers_weight;

            nodes_indices_type m_dfs_indices, m_bfs_indices, m_bfs_levels, m_storage_indices,
                m_any_order_levels;

            basins_type m_basins;
            std::vector<size_type> m_outlets;
            std::vector<size_type> m_pits;

            std::unordered_set<size_type> m_base_levels;
            xt::xarray<bool> m_mask;
            bool m_mask_initialized = false;

            template <class FG, class FR>
            friend class flow_router_impl;

            template <class FG, class SR>
            friend class sink_resolver_impl;

            template <class FG, class OP, class Tag>
            friend class flow_operator_impl;
        };

        template <class G, class S>
        void flow_graph_impl<G, S, flow_graph_fixed_array_tag>::compute_donors()
        {
            m_donors_count.fill(0);

            for (size_type i = 0; i < size(); ++i)
            {
                if (m_receivers(i, 0) != i)
                {
                    // TODO: make it work with multiple flow
                    auto irec = m_receivers(i, 0);
                    m_donors(irec, m_donors_count(irec)++) = i;
                }
            }
        }

        /*
         * Perform depth-first search in the downstream->upstream direction
         * and store the node indices for faster graph traversal.
         */
        template <class G, class S>
        void flow_graph_impl<G, S, flow_graph_fixed_array_tag>::compute_dfs_indices_bottomup()
        {
            size_type nstack = 0;

            std::stack<size_type> tmp;

            for (size_type i = 0; i < size(); ++i)
            {
                if (m_receivers(i, 0) == i)
                {
                    tmp.push(i);
                    m_dfs_indices(nstack++) = i;
                }

                while (!tmp.empty())
                {
                    size_type istack = tmp.top();
                    tmp.pop();

                    for (size_type k = 0; k < m_donors_count(istack); ++k)
                    {
                        const auto idonor = m_donors(istack, k);
                        if (idonor != istack)
                        {
                            m_dfs_indices(nstack++) = idonor;
                            tmp.push(idonor);
                        }
                    }
                }
            }

            assert(nstack == size());
        }

        /*
         * Perform breadth-first search on the flow graph.
         */
        template <class G, class S>
        void flow_graph_impl<G, S, flow_graph_fixed_array_tag>::compute_bfs_indices_bottomup()
        {
            std::vector<std::uint8_t> visited(m_grid.size(), std::uint8_t(0));
            std::vector<size_type> levels(m_grid.size() + 1, 0);
            size_type nstack = 0;
            size_type level = 0;
            bool skip;

            for (size_type i = 0; i < size(); ++i)
                if (m_receivers(i, 0) == i)
                {
                    m_bfs_indices(nstack++) = i;
                }

            levels[level++] = 0;
            levels[level++] = nstack;

            while (nstack < size())
            {
                for (size_type i = levels[level - 2]; i < levels[level - 1]; ++i)
                {
                    auto node_idx = m_bfs_indices(i);
                    visited[node_idx] = 1;

                    for (size_type k = 0; k < m_donors_count(node_idx); ++k)
                    {
                        auto donor_idx = m_donors(node_idx, k);

                        skip = visited[donor_idx] > 0;
                        if (skip)
                            continue;

                        for (std::size_t rcv_idx = 0; rcv_idx < m_receivers_count(donor_idx);
                             ++rcv_idx)
                        {
                            if (visited[m_receivers(donor_idx, rcv_idx)] != 1)
                            {
                                skip = true;
                                break;
                            }
                        }

                        if (!skip && visited[donor_idx] == 0)
                        {
                            m_bfs_indices(nstack++) = donor_idx;
                            visited[donor_idx] = 2;
                        }
                    }
                }

                for (size_type i = levels[level - 2]; i < levels[level - 1]; ++i)
                    visited[m_bfs_indices(i)] = 1;

                levels[level++] = nstack;
            }

            m_bfs_levels = xt::adapt(levels, { level });
        }

        /*
         * Perform depth-first search in the upstream->dowstream direction
         * and store the node indices for faster graph traversal.
         *
         * Note: indices are stored in the downstream->upstream direction!
         */
        template <class G, class S>
        void flow_graph_impl<G, S, flow_graph_fixed_array_tag>::compute_dfs_indices_topdown()
        {
            size_type nstack = 0;
            std::stack<size_type> tmp;

            std::vector<size_type> visited_count(size(), 0);

            for (size_type i = 0; i < size(); ++i)
            {
                if (m_donors_count(i) == 0)
                {
                    tmp.push(i);
                }

                while (!tmp.empty())
                {
                    size_type istack = tmp.top();
                    tmp.pop();

                    m_dfs_indices(nstack++) = istack;

                    for (size_type k = 0; k < m_receivers_count(istack); ++k)
                    {
                        const auto irec = m_receivers(istack, k);
                        visited_count[irec]++;

                        if (visited_count[irec] == m_donors_count(irec))
                        {
                            tmp.push(irec);
                        }
                    }
                }
            }

            assert(nstack == size());

            // downstream->upstream order
            std::reverse(m_dfs_indices.begin(), m_dfs_indices.end());
        }

        template <class G, class S>
        void flow_graph_impl<G, S, flow_graph_fixed_array_tag>::compute_basins()
        {
            size_type current_basin = static_cast<size_type>(-1);
            size_type no_basin = std::numeric_limits<size_type>::max();

            m_outlets.clear();

            for (const auto& inode : nodes_indices_bottomup())
            {
                if (is_masked(inode))
                {
                    m_basins(inode) = no_basin;
                    continue;
                }

                // outlet node has only one receiver: itself
                if (inode == m_receivers(inode, 0))
                {
                    m_outlets.push_back(inode);
                    current_basin++;
                }

                m_basins(inode) = current_basin;
            }

            assert(m_outlets.size() == current_basin + 1);
        }

        template <class G, class S>
        auto flow_graph_impl<G, S, flow_graph_fixed_array_tag>::pits()
            -> const std::vector<size_type>&
        {
            m_pits.clear();

            for (const auto outlet : m_outlets)
            {
                if (!is_base_level(outlet))
                {
                    m_pits.push_back(outlet);
                }
            }

            return m_pits;
        }

        template <class G, class S>
        template <class T>
        void flow_graph_impl<G, S, flow_graph_fixed_array_tag>::accumulate(data_array_type& acc,
                                                                           T&& src) const
        {
            // TODO: assert(acc.shape() == m_grid.shape()) with compatible shape types
            auto src_arr = xt::broadcast(std::forward<T>(src), m_grid.shape());

            // re-init accumulated values
            acc.fill(0);

            auto nodes_indices = nodes_indices_bottomup();

            for (auto inode_ptr = nodes_indices.rbegin(); inode_ptr != nodes_indices.rend();
                 ++inode_ptr)
            {
                const auto inode = *inode_ptr;

                acc.flat(inode) += m_grid.nodes_areas(inode) * src_arr(inode);

                for (size_type r = 0; r < m_receivers_count[inode]; ++r)
                {
                    size_type ireceiver = m_receivers(inode, r);
                    if (ireceiver != inode)
                    {
                        acc.flat(ireceiver) += acc.flat(inode) * m_receivers_weight(inode, r);
                    }
                }
            }
        }

        template <class G, class S>
        template <class T>
        auto flow_graph_impl<G, S, flow_graph_fixed_array_tag>::accumulate(T&& src) const
            -> data_array_type
        {
            data_array_type acc = data_array_type::from_shape(m_grid.shape());
            accumulate(acc, std::forward<T>(src));
            return acc;
        }
    }
}


#endif  // FASTSCAPELIB_FLOW_GRAPH_IMPL_H_
