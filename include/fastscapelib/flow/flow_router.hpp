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

#include <stack>

#include "fastscapelib/algo/flow_routing.hpp"
#include "fastscapelib/flow/flow_graph_impl.hpp"
#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{

    namespace detail
    {

        /**
         * Common implementation for all flow routing methods.
         *
         * @tparam FG The flow graph implementation type.
         * @tparam FR The flow router type.
         */
        template <class FG, class FR>
        class flow_router_impl_base
        {
        public:
            using graph_impl_type = FG;
            using router_type = FR;

        protected:
            flow_router_impl_base(graph_impl_type& graph_impl, const router_type& router)
                : m_graph_impl(graph_impl)
                , m_router(router){};

            ~flow_router_impl_base() = default;

            graph_impl_type& m_graph_impl;
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
         * @tparam FG The flow graph implememtation type.
         * @tparam FR The flow router type.
         */
        template <class FG, class FR>
        class flow_router_impl : public flow_router_impl_base<FG, FR>
        {
        public:
            using graph_impl_type = FG;
            using flow_router_type = FR;
            using base_type = flow_router_impl_base<graph_impl_type, flow_router_type>;

            using data_array_type = typename graph_impl_type::data_array_type;

            // we don't want to instantiate a generic implementation
            // -> only support calling a specialized class template constructor
            flow_router_impl(graph_impl_type& graph_impl, const flow_router_type& router) = delete;

            void route1(const data_array_type& /*elevation*/){};
            void route2(const data_array_type& /*elevation*/){};
        };
    }


    struct single_flow_router
    {
        using flow_graph_impl_tag = detail::flow_graph_fixed_array_tag;
        static constexpr bool is_single = true;
    };


    namespace detail
    {

        template <class FG>
        class flow_router_impl<FG, single_flow_router>
            : public flow_router_impl_base<FG, single_flow_router>
        {
        public:
            using graph_impl_type = FG;
            using base_type = flow_router_impl_base<graph_impl_type, single_flow_router>;

            using data_array_type = typename graph_impl_type::data_array_type;

            flow_router_impl(graph_impl_type& graph_impl, const single_flow_router& router)
                : base_type(graph_impl, router)
            {
                // single flow -> constant weights and receiver_count
                this->m_graph_impl.m_receivers_count.fill(1);
                auto weights = xt::col(this->m_graph_impl.m_receivers_weight, 0);
                weights.fill(1.);
            };

            void route1(const data_array_type& elevation)
            {
                using neighbors_type = typename graph_impl_type::grid_type::neighbors_type;

                double slope, slope_max;
                neighbors_type neighbors;

                auto& grid = this->m_graph_impl.grid();
                auto& donors = this->m_graph_impl.m_donors;
                auto& donors_count = this->m_graph_impl.m_donors_count;
                auto& receivers = this->m_graph_impl.m_receivers;
                auto& dist2receivers = this->m_graph_impl.m_receivers_distance;

                donors_count.fill(0);

                for (auto i : grid.nodes_indices())
                {
                    receivers(i, 0) = i;
                    dist2receivers(i, 0) = 0;
                    slope_max = std::numeric_limits<double>::min();

                    if (grid.status_at_nodes().data()[i] == node_status::fixed_value_boundary)
                    {
                        continue;
                    }

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

                    // fastpath for single flow
                    auto irec = receivers(i, 0);
                    donors(irec, donors_count(irec)++) = i;
                }

                this->m_graph_impl.compute_dfs_indices_downup();
            };

            void route2(const data_array_type& /*elevation*/){};
        };
    }


    struct multiple_flow_router
    {
        using flow_graph_impl_tag = detail::flow_graph_fixed_array_tag;
        static constexpr bool is_single = false;
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
            using graph_impl_type = FG;
            using base_type = flow_router_impl_base<graph_impl_type, multiple_flow_router>;

            using data_array_type = typename graph_impl_type::data_array_type;

            static constexpr size_t n_receivers = graph_impl_type::grid_type::n_neighbors_max();

            flow_router_impl(graph_impl_type& graph_impl, const multiple_flow_router& router)
                : base_type(graph_impl, router){};

            void route1(const data_array_type& /*elevation*/){};
            void route2(const data_array_type& /*elevation*/){};
        };
    }
}

#endif
