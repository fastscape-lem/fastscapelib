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
            using index_type = typename graph_type::index_type;

            using donors_type = typename graph_type::donors_type;
            using donors_count_type = typename graph_type::donors_count_type;

            using receivers_type = typename graph_type::receivers_type;
            using receivers_count_type = typename graph_type::receivers_count_type;
            using receivers_distance_type = typename graph_type::receivers_distance_type;
            using receivers_weight_type = typename graph_type::receivers_weight_type;

            using stack_type = typename graph_type::stack_type;

            flow_router_impl_base(graph_type& graph, router_type& router)
                : m_graph(graph)
                , m_router(router){};

            ~flow_router_impl_base() = default;

            donors_type& donors()
            {
                return m_graph.m_donors;
            };
            donors_count_type& donors_count()
            {
                return m_graph.m_donors_count;
            };

            receivers_type& receivers()
            {
                return m_graph.m_receivers;
            };
            receivers_count_type& receivers_count()
            {
                return m_graph.m_receivers_count;
            };
            receivers_distance_type& receivers_distance()
            {
                return m_graph.m_receivers_distance;
            };
            receivers_weight_type& receivers_weight()
            {
                return m_graph.m_receivers_weight;
            };

            stack_type& dfs_stack()
            {
                return m_graph.m_dfs_stack;
            };

            graph_type& m_graph;
            router_type& m_router;
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

            using elevation_type = typename graph_type::elevation_type;

            static constexpr size_t n_receivers = 0;

            flow_router_impl(graph_type& graph, flow_router_type& router) = delete;

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

            flow_router_impl(graph_type& graph, single_flow_router& router)
                : base_type(graph, router){};

            void route1(const elevation_type& elevation)
            {
                using neighbors_type = typename graph_type::grid_type::neighbors_type;

                double slope, slope_max;
                neighbors_type neighbors;

                auto& grid = this->m_graph.grid();
                auto& donors = this->donors();
                auto& donors_count = this->donors_count();
                auto& receivers = this->receivers();
                auto& dist2receivers = this->receivers_distance();

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

                this->receivers_count().fill(1);

                auto weights = xt::col(this->receivers_weight(), 0);
                weights.fill(1.);

                compute_dfs_stack();
            };

            void route2(const elevation_type& /*elevation*/){};

        private:
            using index_type = typename base_type::index_type;
            using stack_type = typename base_type::stack_type;
            using donors_count_type = typename base_type::donors_count_type;
            using donors_type = typename base_type::donors_type;

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
                const auto& receivers = this->receivers();
                const auto& donors = this->donors();
                const auto& donors_count = this->donors_count();

                auto& stack = this->dfs_stack();
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
            using base_type = flow_router_impl_base<graph_type, single_flow_router>;

            using elevation_type = typename graph_type::elevation_type;

            static constexpr size_t n_receivers = graph_type::grid_type::max_neighbors();

            flow_router_impl(graph_type& graph, single_flow_router& router)
                : base_type(graph, router){};

            void route1(const elevation_type& /*elevation*/){};
            void route2(const elevation_type& /*elevation*/){};
        };
    }
}

#endif
