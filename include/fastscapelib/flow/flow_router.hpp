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

#include <cmath>
#include <map>
#include <stack>
#include <string>

#include "xtensor/xview.hpp"

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

            using impl_embedded_type = std::map<std::string, graph_impl_type>;

            const impl_embedded_type& impl_embedded() const
            {
                return m_impl_embedded;
            }

        protected:
            flow_router_impl_base(graph_impl_type& graph_impl, const router_type& router)
                : m_graph_impl(graph_impl)
                , m_router(router){};

            ~flow_router_impl_base() = default;

            graph_impl_type& m_graph_impl;
            const router_type& m_router;
            impl_embedded_type m_impl_embedded;
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
                        slope = (elevation.flat(i) - elevation.flat(n.idx)) / n.distance;

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

                this->m_graph_impl.compute_dfs_indices_bottomup();
            };

            void route2(const data_array_type& /*elevation*/){};
        };
    }


    struct multi_flow_router
    {
        using flow_graph_impl_tag = detail::flow_graph_fixed_array_tag;
        static constexpr bool is_single = false;
        double slope_exp = 1.0;
    };


    namespace detail
    {

        template <class FG>
        class flow_router_impl<FG, multi_flow_router>
            : public flow_router_impl_base<FG, multi_flow_router>
        {
        public:
            using graph_impl_type = FG;
            using base_type = flow_router_impl_base<graph_impl_type, multi_flow_router>;

            using data_array_type = typename graph_impl_type::data_array_type;

            static constexpr size_t n_receivers = graph_impl_type::grid_type::n_neighbors_max();

            flow_router_impl(graph_impl_type& graph_impl, const multi_flow_router& router)
                : base_type(graph_impl, router){};

            void route1(const data_array_type& elevation)
            {
                using neighbors_type = typename graph_impl_type::grid_type::neighbors_type;
                using nrec_type = typename graph_impl_type::grid_type::neighbors_count_type;

                double slope;
                double weight, weights_sum;
                neighbors_type neighbors;
                nrec_type nrec;

                auto& grid = this->m_graph_impl.grid();
                auto& donors = this->m_graph_impl.m_donors;
                auto& donors_count = this->m_graph_impl.m_donors_count;
                auto& receivers = this->m_graph_impl.m_receivers;
                auto& receivers_count = this->m_graph_impl.m_receivers_count;
                auto& receivers_weight = this->m_graph_impl.m_receivers_weight;
                auto& dist2receivers = this->m_graph_impl.m_receivers_distance;

                donors_count.fill(0);

                for (auto i : grid.nodes_indices())
                {
                    if (grid.status_at_nodes().data()[i] == node_status::fixed_value_boundary)
                    {
                        receivers_count(i) = 1;
                        receivers(i, 0) = i;
                        receivers_weight(i, 0) = 0;
                        dist2receivers(i, 0) = 0;
                        continue;
                    }

                    nrec = 0;
                    weights_sum = 0;

                    for (auto n : grid.neighbors(i, neighbors))
                    {
                        if (elevation.flat(i) > elevation.flat(n.idx))
                        {
                            slope = (elevation.flat(i) - elevation.flat(n.idx)) / n.distance;

                            receivers(i, nrec) = n.idx;
                            dist2receivers(i, nrec) = n.distance;

                            weight = std::pow(slope, this->m_router.slope_exp);
                            weights_sum += weight;
                            receivers_weight(i, nrec) = weight;

                            // update donors (note: not thread safe if later parallelization)
                            donors(n.idx, donors_count(n.idx)++) = i;

                            nrec++;
                        }
                    }

                    receivers_count(i) = nrec;

                    // normalize weights
                    for (auto j = 0; j < nrec; j++)
                    {
                        receivers_weight(i, j) /= weights_sum;
                    }
                }

                // DFS upstream->downstream so that it works with multi-directions flow
                this->m_graph_impl.compute_dfs_indices_topdown();
            };

            void route2(const data_array_type& /*elevation*/){};
        };
    }


    class singlemulti_flow_router
    {
    public:
        using flow_graph_impl_tag = detail::flow_graph_fixed_array_tag;
        static constexpr bool is_single = false;
        double slope_exp = 1.0;
        bool store_single_unresolved = false;
        bool store_single_resolved = false;
    };


    namespace detail
    {
        template <class FG>
        class flow_router_impl<FG, singlemulti_flow_router>
            : public flow_router_impl_base<FG, singlemulti_flow_router>
        {
        public:
            using graph_impl_type = FG;
            using base_type = flow_router_impl_base<graph_impl_type, singlemulti_flow_router>;

            using data_array_type = typename graph_impl_type::data_array_type;

            flow_router_impl(graph_impl_type& graph_impl, const singlemulti_flow_router& router)
                : base_type(graph_impl, router)
                , m_single_router_impl(graph_impl, single_flow_router())
                , m_multi_router_impl(graph_impl, { router.slope_exp })
                , m_store_single_unresolved(router.store_single_unresolved)
                , m_store_single_resolved(router.store_single_resolved)
            {
                auto& grid = graph_impl.grid();
                single_flow_router srouter;

                if (m_store_single_unresolved)
                {
                    this->m_impl_embedded.insert(
                        { "single_unresolved", graph_impl_type(grid, srouter) });
                }
                if (m_store_single_resolved)
                {
                    this->m_impl_embedded.insert(
                        { "single_resolved", graph_impl_type(grid, srouter) });
                }
            };

            void route1(const data_array_type& elevation)
            {
                // weights and receiver_count not reset in `single_flow_router::route1`
                // -> do it here
                this->m_graph_impl.m_receivers_count.fill(1);
                auto weights = xt::col(this->m_graph_impl.m_receivers_weight, 0);
                weights = 1;

                m_single_router_impl.route1(elevation);

                if (m_store_single_unresolved)
                {
                    copy_graph_impl("single_unresolved");
                }
            };

            void route2(const data_array_type& elevation)
            {
                if (m_store_single_resolved)
                {
                    copy_graph_impl("single_resolved");
                }

                m_multi_router_impl.route1(elevation);
            };

        private:
            bool m_store_single_unresolved;
            bool m_store_single_resolved;

            flow_router_impl<FG, single_flow_router> m_single_router_impl;
            flow_router_impl<FG, multi_flow_router> m_multi_router_impl;

            /*
             * Copy (multiple flow) graph implementation data to one of the
             * embedded (single flow) graphs
             */
            void copy_graph_impl(std::string key)
            {
                const auto& mgraph = this->m_graph_impl;
                auto& sgraph = this->m_impl_embedded.at(key);

                sgraph.m_receivers_count = mgraph.m_receivers_count;

                auto receivers_col = xt::col(sgraph.m_receivers, 0);
                receivers_col = xt::col(mgraph.m_receivers, 0);
                auto receivers_distance_col = xt::col(sgraph.m_receivers_distance, 0);
                receivers_distance_col = xt::col(mgraph.m_receivers_distance, 0);
                auto receivers_weight_col = xt::col(sgraph.m_receivers_weight, 0);
                receivers_weight_col = xt::col(mgraph.m_receivers_weight, 0);

                sgraph.m_donors_count = mgraph.m_donors_count;

                auto donors_col = xt::col(sgraph.m_donors, 0);
                donors_col = xt::col(mgraph.m_donors, 0);

                sgraph.m_dfs_indices = mgraph.m_dfs_indices;
            }
        };
    }
}

#endif
