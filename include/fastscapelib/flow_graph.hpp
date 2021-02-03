/**
 * @brief Class used to compute or follow flow routes on
 * the topographic surface.
 * 
 * It provides a single API to chain flow_router and
 * sink_resolver methods and access their results.
 *
 */
#ifndef FASTSCAPELIB_FLOW_GRAPH_H
#define FASTSCAPELIB_FLOW_GRAPH_H

#include "fastscapelib/flow_router.hpp"
#include "fastscapelib/sink_resolver.hpp"
#include "fastscapelib/xtensor_utils.hpp"

#include <stdexcept>


namespace fastscapelib
{

    template <class FG>
    class sink_resolver;

    /**
     * Main class used to compute or follow flow routes on
     * the topographic surface.
     *
     * @tparam G The grid type.
     * @tparam elev_t Type used to store elevation values.
     * @tparam G The xtensor selector type.
     */
    template <class G, class elev_t, class S = typename G::xt_selector>
    class flow_graph
    {
    public:
        using self_type = flow_graph<G, elev_t, S>;
        using grid_type = G;

        using index_type = typename grid_type::size_type;
        using distance_type = typename grid_type::distance_type;
        using neighbors_count_type = typename grid_type::neighbors_count_type;

        using elevation_type = xt_container_t<S, elev_t, G::xt_ndims>;

        using donors_type = xt_container_t<S, index_type, 2>;
        using donors_count_type = xt_container_t<S, neighbors_count_type, 1>;

        using receivers_type = donors_type;
        using receivers_count_type = donors_count_type;
        using receivers_weight_type = xt_container_t<S, double, 2>;
        using receivers_distance_type = xt_container_t<S, distance_type, 2>;

        using stack_type = xt_container_t<S, index_type, 1>;

        using const_dfs_iterator = const index_type*;
        using const_reverse_dfs_iterator = std::reverse_iterator<const index_type*>;

        using flow_router_ptr = std::unique_ptr<flow_router<self_type>>;
        using sink_resolver_ptr = std::unique_ptr<sink_resolver<self_type>>;

        flow_graph(G& grid, flow_router_ptr router, sink_resolver_ptr resolver)
            : m_grid(grid),
              p_flow_router(std::move(router)),
              p_sink_resolver(std::move(resolver))
        {
            shape_type receivers_shape = {grid.size(), grid_type::max_neighbors()};
            shape_type donors_shape = {grid.size(), grid_type::max_neighbors()+1};

            m_receivers = xt::ones<index_type>(receivers_shape) * -1;
            m_receivers_count = xt::zeros<index_type>({grid.size()});
            m_receivers_distance = xt::ones<distance_type>(receivers_shape) * -1;
            m_receivers_weight = xt::zeros<double>(receivers_shape);

            m_donors = xt::ones<index_type>(donors_shape) * -1;
            m_donors_count = xt::zeros<index_type>({grid.size()});

            m_dfs_stack = xt::ones<index_type>({grid.size()}) * -1;
        }

        const elevation_type& update_routes(const elevation_type& elevation)
        {
            const auto& modified_elevation = p_sink_resolver->resolve_before_route(elevation, *this);
            p_flow_router->route(modified_elevation, *this);
            const auto& final_elevation = p_sink_resolver->resolve_after_route(modified_elevation, *this);

            return final_elevation;
        }

        G& grid() { return m_grid; };

        index_type size() { return m_grid.size(); };

        const receivers_type& receivers() { return m_receivers; };

        const receivers_count_type& receivers_count() { return m_receivers_count; };

        const receivers_distance_type& receivers_distance() { return m_receivers_distance; };

        const receivers_weight_type& receivers_weight() { return m_receivers_weight; };

        const donors_type& donors() { return m_donors; };

        const donors_count_type& donors_count() { return m_donors_count; };

        const stack_type& dfs_stack() { return m_dfs_stack; };

        const_dfs_iterator dfs_cbegin() { return m_dfs_stack.cbegin(); };

        const_dfs_iterator dfs_cend() { return m_dfs_stack.cend(); };

        const_reverse_dfs_iterator dfs_crbegin() { return m_dfs_stack.crbegin(); };

        const_reverse_dfs_iterator dfs_crend() { return m_dfs_stack.crend(); };

    private:
    
        using shape_type = typename xt_container_t<S, index_type, 2>::shape_type;

        G& m_grid;

        donors_type m_donors;
        donors_count_type m_donors_count;
        
        receivers_type m_receivers;
        receivers_count_type m_receivers_count;
        receivers_distance_type m_receivers_distance;
        receivers_weight_type m_receivers_weight;

        stack_type m_dfs_stack;

        flow_router_ptr p_flow_router;
        sink_resolver_ptr p_sink_resolver;

        friend class flow_router<self_type>;
        friend class sink_resolver<self_type>;
    };
}

#endif