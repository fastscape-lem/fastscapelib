/**
 * @brief Classes used to compute flow routes on
 * the topographic surface.
 *
 * ``flow_router`` is meant to be used in combination
 * with possible ``sink_resolver`` through a ``flow_graph``.
 *
 */
#ifndef FASTSCAPELIB_FLOW_FLOW_ROUTER_H
#define FASTSCAPELIB_FLOW_FLOW_ROUTER_H

#include "fastscapelib/algo/flow_routing.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{

    namespace detail
    {

        /**
         * Common implementation for all flow routing methods.
         *
         * @tparam FG The flow graph type.
         * @tparam FR The flow router type.
         */
        template <class FG, class FR>
        class flow_router_impl_base
        {
        public:
            using graph_type = FG;
            using router_type = FR;

        protected:
            flow_router_impl_base(graph_type& graph, const router_type& router)
                : m_graph(graph)
                , m_router(router){};

            ~flow_router_impl_base() = default;

            graph_type& m_graph;
            const router_type& m_router;
        };


        /**
         * Flow routing implementation.
         *
         * This class is used via template specialization (one for each
         * flow routing method).
         *
         * The declaration for the generic case here contains the minimum that
         * should be (re)implemented in specialized template classes.
         *
         * @tparam FG The flow graph type.
         * @tparam FR The flow router type.
         */
        template <class FG, class FR>
        class flow_router_impl : public flow_router_impl_base<FG, FR>
        {
        public:
            using graph_type = FG;
            using flow_router_type = FR;
            using base_type = flow_router_impl_base<graph_type, flow_router_type>;

            using elevation_type = typename graph_type::elevation_type;

            static constexpr size_t n_receivers = 0;

            // we don't want to instantiate a generic implementation
            // -> only support calling a specialized class template constructor
            flow_router_impl(graph_type& graph, const flow_router_type& router) = delete;

            void route1(const elevation_type& /*elevation*/){};
            void route2(const elevation_type& /*elevation*/){};
        };
    }


    struct single_flow_router
    {
    };


    namespace detail
    {

        template <class FG>
        class flow_router_impl<FG, single_flow_router>
            : public flow_router_impl_base<FG, single_flow_router>
        {
        public:
            using graph_type = FG;
            using base_type = flow_router_impl_base<graph_type, single_flow_router>;

            using elevation_type = typename graph_type::elevation_type;

            static constexpr size_t n_receivers = 1;

            flow_router_impl(graph_type& graph, const single_flow_router& router)
                : base_type(graph, router){};

            void route1(const elevation_type& elevation)
            {
                using neighbors_type = typename graph_type::grid_type::neighbors_type;

                double slope, slope_max;
                neighbors_type neighbors;

                auto& grid = this->m_graph.grid();
                auto& donors = this->m_graph.m_donors;
                auto& donors_count = this->m_graph.m_donors_count;
                auto& receivers = this->m_graph.m_receivers;
                auto& dist2receivers = this->m_graph.m_receivers_distance;

                donors_count.fill(0);

                for (auto i : grid.nodes_indices())
                {
                    receivers(i, 0) = i;
                    dist2receivers(i, 0) = 0;
                    slope_max = std::numeric_limits<double>::min();

                    for (auto n : grid.neighbors(i, neighbors))
                    {
                        slope = (elevation.data()[i] - elevation.data()[n.idx]) / n.distance;

                        if (slope > slope_max)
                        {
                            slope_max = slope;
                            receivers(i, 0) = n.idx;
                            dist2receivers(i, 0) = n.distance;
                        }
                    }
                    donors(receivers(i, 0), donors_count(receivers(i, 0))++) = i;
                }

                this->m_graph.m_receivers_count.fill(1);

                auto weights = xt::col(this->m_graph.m_receivers_weight, 0);
                weights.fill(1.);

                compute_dfs_stack();
            };

            void route2(const elevation_type& /*elevation*/){};

        private:
            using index_type = typename graph_type::index_type;
            using stack_type = typename graph_type::stack_type;
            using donors_count_type = typename graph_type::donors_count_type;
            using donors_type = typename graph_type::donors_type;

            void add2stack(index_type& nstack,
                           stack_type& stack,
                           const donors_count_type& ndonors,
                           const donors_type& donors,
                           const index_type inode)
            {
                for (index_type k = 0; k < ndonors(inode); ++k)
                {
                    const auto idonor = donors(inode, k);
                    if (idonor != inode)
                    {
                        stack(nstack++) = idonor;
                        add2stack(nstack, stack, ndonors, donors, idonor);
                    }
                }
            }

            void compute_dfs_stack()
            {
                const auto& receivers = this->m_graph.m_receivers;
                const auto& donors = this->m_graph.m_donors;
                const auto& donors_count = this->m_graph.m_donors_count;

                auto& stack = this->m_graph.m_dfs_stack;
                index_type nstack = 0;

                for (index_type i = 0; i < this->m_graph.size(); ++i)
                {
                    if (receivers(i, 0) == i)
                    {
                        stack(nstack++) = i;
                        add2stack(nstack, stack, donors_count, donors, i);
                    }
                }
            };
        };
    }


    struct multiple_flow_router
    {
        double p1;
        double p2;
    };


    namespace detail
    {

        /**
         * TODO: not yet operational.
         *
         */
        template <class FG>
        class flow_router_impl<FG, multiple_flow_router>
            : public flow_router_impl_base<FG, multiple_flow_router>
        {
        public:
            using graph_type = FG;
            using base_type = flow_router_impl_base<graph_type, multiple_flow_router>;

            using elevation_type = typename graph_type::elevation_type;

            static constexpr size_t n_receivers = graph_type::grid_type::max_neighbors();

            flow_router_impl(graph_type& graph, const multiple_flow_router& router)
                : base_type(graph, router){};

            void route1(const elevation_type& /*elevation*/){};
            void route2(const elevation_type& /*elevation*/){};
        };
    }
}

#endif
