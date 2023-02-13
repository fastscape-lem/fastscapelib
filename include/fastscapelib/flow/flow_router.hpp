#ifndef FASTSCAPELIB_FLOW_FLOW_ROUTER_H
#define FASTSCAPELIB_FLOW_FLOW_ROUTER_H


#include "fastscapelib/flow/flow_graph_impl.hpp"
#include "fastscapelib/flow/flow_operator.hpp"
#include "fastscapelib/grid/base.hpp"


namespace fastscapelib
{

    /**
     * Single direction direction flow router operator.
     *
     * This flow operator routes all the flow passing through a grid node
     * towards its neighbors node of steepest slope.
     */
    class single_flow_router : public flow_operator
    {
    public:
        inline std::string name() const noexcept override
        {
            return "single_flow_router";
        }

        static constexpr bool graph_updated = true;
        static constexpr flow_direction out_flowdir = flow_direction::single;
    };


    namespace detail
    {

        /**
         * Single flow router operator implementation.
         */
        template <class FG>
        class flow_operator_impl<FG, single_flow_router, flow_graph_fixed_array_tag>
            : public flow_operator_impl_base<FG, single_flow_router>
        {
        public:
            using graph_impl_type = FG;
            using base_type = flow_operator_impl_base<FG, single_flow_router>;
            using data_array_type = typename graph_impl_type::data_array_type;

            flow_operator_impl(std::shared_ptr<single_flow_router> ptr)
                : base_type(std::move(ptr)){};

            void apply(graph_impl_type& graph_impl, data_array_type& elevation)
            {
                // single flow optimization
                graph_impl.m_receivers_count.fill(1);
                auto weights = xt::col(graph_impl.m_receivers_weight, 0);
                weights.fill(1.);

                using neighbors_type = typename graph_impl_type::grid_type::neighbors_type;

                double slope, slope_max;
                neighbors_type neighbors;

                auto& grid = graph_impl.grid();
                auto& donors = graph_impl.m_donors;
                auto& donors_count = graph_impl.m_donors_count;
                auto& receivers = graph_impl.m_receivers;
                auto& dist2receivers = graph_impl.m_receivers_distance;

                donors_count.fill(0);

                for (auto i : grid.node_indices())
                {
                    receivers(i, 0) = i;
                    dist2receivers(i, 0) = 0;
                    slope_max = std::numeric_limits<double>::min();

                    if (grid.status_at_nodes().flat(i) == node_status::fixed_value_boundary)
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

                graph_impl.compute_dfs_indices_bottomup();
            }
        };
    }


    /**
     * Multiple direction flow router operator.
     *
     * This flow operator partitions the flow passing through a grid node among
     * its downslope neighbor nodes. Flow partitioning is proportional to the
     * local slope between a node and its neighbors (power relationship with a
     * fixed exponent parameter).
     */
    class multi_flow_router : public flow_operator
    {
    public:
        inline std::string name() const noexcept override
        {
            return "multi_flow_router";
        }

        static constexpr bool graph_updated = true;
        static constexpr flow_direction out_flowdir = flow_direction::multi;

        /**
         * Create a new multi flow router operator.
         *
         * @param slope_exp The flow partition slope exponent.
         */
        multi_flow_router(double slope_exp)
            : m_slope_exp(slope_exp)
        {
        }

        double m_slope_exp = 1.0;
    };


    namespace detail
    {

        /**
         * Multiple direction flow router operator implementation.
         */
        template <class FG>
        class flow_operator_impl<FG, multi_flow_router, flow_graph_fixed_array_tag>
            : public flow_operator_impl_base<FG, multi_flow_router>
        {
        public:
            using graph_impl_type = FG;
            using base_type = flow_operator_impl_base<FG, multi_flow_router>;
            using data_array_type = typename graph_impl_type::data_array_type;

            flow_operator_impl(std::shared_ptr<multi_flow_router> ptr)
                : base_type(std::move(ptr)){};

            void apply(graph_impl_type& graph_impl, data_array_type& elevation)
            {
                using neighbors_type = typename graph_impl_type::grid_type::neighbors_type;
                using nrec_type = typename graph_impl_type::grid_type::neighbors_count_type;

                double slope;
                double weight, weights_sum;
                neighbors_type neighbors;
                nrec_type nrec;

                auto& grid = graph_impl.grid();
                auto& donors = graph_impl.m_donors;
                auto& donors_count = graph_impl.m_donors_count;
                auto& receivers = graph_impl.m_receivers;
                auto& receivers_count = graph_impl.m_receivers_count;
                auto& receivers_weight = graph_impl.m_receivers_weight;
                auto& dist2receivers = graph_impl.m_receivers_distance;

                donors_count.fill(0);

                for (auto i : grid.node_indices())
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

                            weight = std::pow(slope, this->m_op_ptr->m_slope_exp);
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
                graph_impl.compute_dfs_indices_topdown();
            }
        };
    }
}

#endif
