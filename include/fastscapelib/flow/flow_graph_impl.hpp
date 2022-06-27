#ifndef FASTSCAPELIB_FLOW_GRAPH_IMPL_H_
#define FASTSCAPELIB_FLOW_GRAPH_IMPL_H_

#include <array>

#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{
    namespace detail
    {

        template <class G, class S, class Tag>
        class flow_graph_impl;


        struct flow_graph_fixed_array_tag
        {
        };


        template <class G, class S>
        class flow_graph_impl<G, S, flow_graph_fixed_array_tag>
        {
        public:
            using grid_type = G;

            using grid_data_type = typename grid_type::grid_data_type;

            flow_graph_impl(grid_type& grid)
                : m_grid(grid){};

            const xt_array_t<S, grid_data_type> test() const
            {
                using index_type = typename grid_type::size_type;
                using shape_type = std::array<index_type, 2>;
                const shape_type receivers_shape = { m_grid.size(), grid_type::n_neighbors_max() };
                return xt::ones<grid_data_type>(receivers_shape) * -1;
            }

        private:
            grid_type& m_grid;
        };

    }
}


#endif  // FASTSCAPELIB_FLOW_GRAPH_IMPL_H_
