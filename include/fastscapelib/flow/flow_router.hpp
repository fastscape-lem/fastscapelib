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

    /**
     * Base class for the implementation of flow routing
     * methods.
     *
     * All derived classes must implement ``route1`` and ``route2`` methods.
     *
     * @tparam FG The flow_graph class.
     */
    template <class FG>
    class flow_router
    {
    public:
        using elevation_type = typename FG::elevation_type;

        // Entity semantics, i.e., a flow graph uses a flow router via this base class.
        // -> avoid incomplete destruction (e.g., there may be members in inherited classes
        // that need to be destroyed)
        virtual ~flow_router() = default;

        // do not assign or copy flow routers using the base class.
        flow_router(const flow_router&) = delete;
        flow_router(flow_router&&) = delete;
        flow_router& operator=(const flow_router&) = delete;
        flow_router& operator=(flow_router&&) = delete;

        virtual void route1(const elevation_type& elevation, FG& fgraph) = 0;
        virtual void route2(const elevation_type& elevation, FG& fgraph) = 0;

    protected:
        using index_type = typename FG::index_type;

        using donors_type = typename FG::donors_type;
        using donors_count_type = typename FG::donors_count_type;

        using receivers_type = typename FG::receivers_type;
        using receivers_count_type = typename FG::receivers_count_type;
        using receivers_distance_type = typename FG::receivers_distance_type;
        using receivers_weight_type = typename FG::receivers_weight_type;

        using stack_type = typename FG::stack_type;

        flow_router() = default;

        donors_type& donors(FG& fgraph)
        {
            return fgraph.m_donors;
        };
        donors_count_type& donors_count(FG& fgraph)
        {
            return fgraph.m_donors_count;
        };

        receivers_type& receivers(FG& fgraph)
        {
            return fgraph.m_receivers;
        };
        receivers_count_type& receivers_count(FG& fgraph)
        {
            return fgraph.m_receivers_count;
        };
        receivers_distance_type& receivers_distance(FG& fgraph)
        {
            return fgraph.m_receivers_distance;
        };
        receivers_weight_type& receivers_weight(FG& fgraph)
        {
            return fgraph.m_receivers_weight;
        };

        stack_type& dfs_stack(FG& fgraph)
        {
            return fgraph.m_dfs_stack;
        };
    };


    /**
     * A flow_router not doing anything.
     *
     * @tparam FG The flow_graph class.
     */
    template <class FG>
    class dummy_flow_router final : public flow_router<FG>
    {
    public:
        using base_type = flow_router<FG>;
        using elevation_type = typename base_type::elevation_type;

        dummy_flow_router() = default;

        virtual ~dummy_flow_router() = default;

        void route1(const elevation_type& /*elevation*/, FG& /*fgraph*/){};
        void route2(const elevation_type& /*elevation*/, FG& /*fgraph*/){};
    };


    /**
     * A flow_router considering only one receiver per
     * grid node.
     *
     * @tparam FG The flow_graph class.
     */
    template <class FG>
    class single_flow_router final : public flow_router<FG>
    {
    public:
        using base_type = flow_router<FG>;
        using elevation_type = typename base_type::elevation_type;

        single_flow_router() = default;

        virtual ~single_flow_router() = default;

        void route1(const elevation_type& elevation, FG& fgraph)
        {
            using neighbors_type = typename FG::grid_type::neighbors_type;

            double slope, slope_max;
            neighbors_type neighbors;

            auto& grid = fgraph.grid();
            auto& donors = this->donors(fgraph);
            auto& donors_count = this->donors_count(fgraph);
            auto& receivers = this->receivers(fgraph);
            auto& dist2receivers = this->receivers_distance(fgraph);

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

            this->receivers_count(fgraph).fill(1);

            auto weights = xt::col(this->receivers_weight(fgraph), 0);
            weights.fill(1.);

            compute_dfs_stack(fgraph);
        };

        void route2(const elevation_type& /*elevation*/, FG& /*fgraph*/){};

    private:
        using index_type = typename flow_router<FG>::index_type;
        using stack_type = typename flow_router<FG>::stack_type;
        using donors_count_type = typename flow_router<FG>::donors_count_type;
        using donors_type = typename flow_router<FG>::donors_type;

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

        void compute_dfs_stack(FG& fgraph)
        {
            const auto& receivers = this->receivers(fgraph);
            const auto& donors = this->donors(fgraph);
            const auto& donors_count = this->donors_count(fgraph);

            auto& stack = this->dfs_stack(fgraph);
            index_type nstack = 0;

            for (index_type i = 0; i < fgraph.size(); ++i)
            {
                if (receivers(i, 0) == i)
                {
                    stack(nstack++) = i;
                    add2stack(nstack, stack, donors_count, donors, i);
                }
            }
        };
    };


    /**
     * A flow_router considering multiple receivers per
     * grid node.
     *
     * @tparam FG The flow_graph class.
     */
    template <class FG>
    class multiple_flow_router final : public flow_router<FG>
    {
    public:
        using base_type = flow_router<FG>;
        using elevation_type = typename base_type::elevation_type;

        multiple_flow_router(double param1, double param2)
            : p1(param1)
            , p2(param2)
        {
        }

        virtual ~multiple_flow_router() = default;

        void route1(const elevation_type& /*elevation*/, FG& /*fgraph*/){};
        void route2(const elevation_type& /*elevation*/, FG& /*fgraph*/){};

    private:
        double p1 = 0., p2 = 0.;
    };
}

#endif
