#ifndef FASTSCAPELIB_FLOW_GRAPH_IMPL_H_
#define FASTSCAPELIB_FLOW_GRAPH_IMPL_H_

#include <array>

#include "xtensor/xbroadcast.hpp"
#include "xtensor/xstrided_view.hpp"
#include "xtensor/xview.hpp"

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
            using xt_selector = S;

            using size_type = typename grid_type::size_type;
            using neighbors_count_type = typename grid_type::neighbors_count_type;

            using data_type = typename grid_type::grid_data_type;
            using data_array_type = xt_array_t<xt_selector, data_type>;

            using donors_type = xt_tensor_t<xt_selector, size_type, 2>;
            using donors_count_type = xt_tensor_t<xt_selector, neighbors_count_type, 1>;
            using receivers_type = donors_type;
            using receivers_count_type = donors_count_type;
            using receivers_distance_type = xt_tensor_t<xt_selector, data_type, 2>;
            using receivers_weight_type = xt_tensor_t<xt_selector, data_type, 2>;
            using dfs_indices_type = xt_tensor_t<S, size_type, 1>;

            // using const_dfs_iterator = const size_type*;
            // using const_reverse_dfs_iterator = std::reverse_iterator<const size_type*>;

            flow_graph_impl(grid_type& grid)
                : m_grid(grid)
            {
                using shape_type = std::array<size_type, 2>;
                const shape_type receivers_shape = { grid.size(), grid_type::n_neighbors_max() };
                const shape_type donors_shape = { grid.size(), grid_type::n_neighbors_max() + 1 };

                m_receivers = xt::ones<size_type>(receivers_shape) * -1;
                m_receivers_count = xt::zeros<size_type>({ grid.size() });
                m_receivers_distance = xt::ones<data_type>(receivers_shape) * -1;
                m_receivers_weight = xt::zeros<data_type>(receivers_shape);

                m_donors = xt::ones<size_type>(donors_shape) * -1;
                m_donors_count = xt::zeros<size_type>({ grid.size() });

                m_dfs_indices = xt::ones<size_type>({ grid.size() }) * -1;
            };

            G& grid()
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

            const dfs_indices_type& dfs_indices() const
            {
                return m_dfs_indices;
            };

            // const_dfs_iterator dfs_cbegin()
            // {
            //     return m_dfs_indices.cbegin();
            // };

            // const_dfs_iterator dfs_cend()
            // {
            //     return m_dfs_indices.cend();
            // };

            // const_reverse_dfs_iterator dfs_crbegin()
            // {
            //     return m_dfs_indices.crbegin();
            // };

            // const_reverse_dfs_iterator dfs_crend()
            // {
            //     return m_dfs_indices.crend();
            // };

            template <class T>
            void accumulate(data_array_type& acc, T&& src) const;

            template <class T>
            data_array_type accumulate(T&& src) const;

        private:
            grid_type& m_grid;

            donors_type m_donors;
            donors_count_type m_donors_count;

            receivers_type m_receivers;
            receivers_count_type m_receivers_count;
            receivers_distance_type m_receivers_distance;
            receivers_weight_type m_receivers_weight;

            dfs_indices_type m_dfs_indices;

            template <class FG, class FR>
            friend class flow_router_impl;

            template <class FG, class SR>
            friend class sink_resolver_impl;
        };

        template <class G, class S>
        template <class T>
        void flow_graph_impl<G, S, flow_graph_fixed_array_tag>::accumulate(data_array_type& acc,
                                                                           T&& src) const
        {
            // TODO: assert(acc.shape() == m_grid.shape()) with compatible shape types
            auto src_arr = xt::broadcast(std::forward<T>(src), m_grid.shape());

            // TODO: safer to flatten acc and src? (case of raster_grid -> 2d arrays)
            // flatten seems to greatly slow down the execution
            // alternative? (check if it would work with different layout, e.g., column major)
            // maybe xt::ravel? -> specify row::major explicitly (maybe faster?)
            // or acc.flat(idx) ? -> check perfs (not available for xbroadcast expression)
            // currently: just use operator() works even for 2d arrays. not sure why?

            // re-init accumulated values
            acc.fill(0);

            for (auto inode = m_dfs_indices.crbegin(); inode != m_dfs_indices.crend(); ++inode)
            {
                acc.flat(*inode) += m_grid.node_area(*inode) * src_arr(*inode);

                for (size_type r = 0; r < m_receivers_count[*inode]; ++r)
                {
                    size_type ireceiver = m_receivers(*inode, r);
                    if (ireceiver != *inode)
                    {
                        acc.flat(ireceiver) += acc.flat(*inode) * m_receivers_weight(*inode, r);
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
