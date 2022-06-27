#ifndef FASTSCAPELIB_FLOW_GRAPH_IMPL_H_
#define FASTSCAPELIB_FLOW_GRAPH_IMPL_H_

#include <array>

#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{
    namespace detail
    {

        template <class FG, class FR>
        class flow_router_impl;


        template <class FG, class SR>
        class sink_resolver_impl;


        template <class G, class S, class Tag>
        class flow_graph_impl;


        struct flow_graph_fixed_array_tag
        {
        };


        template <class G, class S>
        class flow_graph_impl<G, S, flow_graph_fixed_array_tag>
        {
        public:
            using self_type = flow_graph_impl<G, S, flow_graph_fixed_array_tag>;
            using grid_type = G;

            using index_type = typename grid_type::size_type;
            using grid_data_type = typename grid_type::grid_data_type;
            using neighbors_count_type = typename grid_type::neighbors_count_type;

            template <class T>
            using data_type = xt_array_t<S, T>;

            // referenced in flow router and sink resolver implementations
            using elevation_type = xt_array_t<S, typename grid_type::grid_data_type>;

            using donors_type = xt_tensor_t<S, index_type, 2>;
            using donors_count_type = xt_tensor_t<S, neighbors_count_type, 1>;

            using receivers_type = donors_type;
            using receivers_count_type = donors_count_type;
            using receivers_weight_type = xt_tensor_t<S, double, 2>;
            using receivers_distance_type = xt_tensor_t<S, grid_data_type, 2>;

            using stack_type = xt_tensor_t<S, index_type, 1>;

            // using const_dfs_iterator = const index_type*;
            // using const_reverse_dfs_iterator = std::reverse_iterator<const index_type*>;

            flow_graph_impl(grid_type& grid)
                : m_grid(grid)
            {
                using shape_type = std::array<index_type, 2>;
                const shape_type receivers_shape = { grid.size(), grid_type::n_neighbors_max() };
                const shape_type donors_shape = { grid.size(), grid_type::n_neighbors_max() + 1 };

                m_receivers = xt::ones<index_type>(receivers_shape) * -1;
                m_receivers_count = xt::zeros<index_type>({ grid.size() });
                m_receivers_distance = xt::ones<grid_data_type>(receivers_shape) * -1;
                m_receivers_weight = xt::zeros<double>(receivers_shape);

                m_donors = xt::ones<index_type>(donors_shape) * -1;
                m_donors_count = xt::zeros<index_type>({ grid.size() });

                m_dfs_stack = xt::ones<index_type>({ grid.size() }) * -1;
            };

            G& grid()
            {
                return m_grid;
            };

            index_type size() const
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

            const stack_type& dfs_stack() const
            {
                return m_dfs_stack;
            };

            // const_dfs_iterator dfs_cbegin()
            // {
            //     return m_dfs_stack.cbegin();
            // };

            // const_dfs_iterator dfs_cend()
            // {
            //     return m_dfs_stack.cend();
            // };

            // const_reverse_dfs_iterator dfs_crbegin()
            // {
            //     return m_dfs_stack.crbegin();
            // };

            // const_reverse_dfs_iterator dfs_crend()
            // {
            //     return m_dfs_stack.crend();
            // };

            template <class T>
            T accumulate(const T& data) const;

            data_type<double> accumulate(const double& data) const;

        private:
            grid_type& m_grid;

            donors_type m_donors;
            donors_count_type m_donors_count;

            receivers_type m_receivers;
            receivers_count_type m_receivers_count;
            receivers_distance_type m_receivers_distance;
            receivers_weight_type m_receivers_weight;

            stack_type m_dfs_stack;

            template <class FG, class FR>
            friend class flow_router_impl;

            template <class FG, class SR>
            friend class sink_resolver_impl;
        };

        template <class G, class S>
        template <class T>
        auto flow_graph_impl<G, S, flow_graph_fixed_array_tag>::accumulate(const T& data) const -> T
        {
            T acc = xt::zeros_like(data);

            for (auto inode = m_dfs_stack.crbegin(); inode != m_dfs_stack.crend(); ++inode)
            {
                acc(*inode) += m_grid.node_area(*inode) * data.data()[*inode];

                for (index_type r = 0; r < m_receivers_count[*inode]; ++r)
                {
                    index_type ireceiver = m_receivers(*inode, r);
                    if (ireceiver != *inode)
                    {
                        acc(ireceiver) += acc(*inode) * m_receivers_weight(*inode, r);
                    }
                }
            }

            return acc;
        }

        template <class G, class S>
        auto flow_graph_impl<G, S, flow_graph_fixed_array_tag>::accumulate(const double& data) const
            -> data_type<double>
        {
            data_type<double> tmp = xt::ones<double>(m_grid.shape()) * data;
            return accumulate(tmp);
        }
    }
}


#endif  // FASTSCAPELIB_FLOW_GRAPH_IMPL_H_
