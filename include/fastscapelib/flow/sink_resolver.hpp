/**
 * @brief Class used to resolve local depressions of
 * the topographic surface.
 *
 * ``sink_resolver`` is meant to be used in combination
 * with ``flow_router`` through a ``flow_graph``.
 *
 */
#ifndef FASTSCAPELIB_FLOW_SINK_RESOLVER_H
#define FASTSCAPELIB_FLOW_SINK_RESOLVER_H

#include <limits>

#include "fastscapelib/algo/pflood.hpp"
#include "fastscapelib/flow/flow_graph_impl.hpp"
#include "fastscapelib/flow/basin_graph.hpp"


namespace fastscapelib
{

    namespace detail
    {

        /**
         * Common implementation for all sink resolving methods.
         *
         * @tparam FG The flow graph implementation type.
         * @tparam SR The sink resolver type.
         */
        template <class FG, class SR>
        class sink_resolver_impl_base
        {
        public:
            using graph_impl_type = FG;
            using resolver_type = SR;

            sink_resolver_impl_base(graph_impl_type& graph_impl, const resolver_type& resolver)
                : m_graph_impl(graph_impl)
                , m_resolver(resolver){};

        protected:
            graph_impl_type& m_graph_impl;
            const resolver_type& m_resolver;
        };


        /**
         * Sink resolving implementation.
         *
         * This class is used via template specialization (one for each
         * flow routing method).
         *
         * The declaration for the generic case here contains the minimum that
         * should be (re)implemented in specialized template classes.
         *
         * @tparam FG The flow graph implementation type.
         * @tparam SR The sink resolver type.
         */
        template <class FG, class SR>
        class sink_resolver_impl : public sink_resolver_impl_base<FG, SR>
        {
        public:
            using graph_impl_type = FG;
            using resolver_type = SR;
            using base_type = sink_resolver_impl_base<graph_impl_type, resolver_type>;

            using data_array_type = typename graph_impl_type::data_array_type;

            sink_resolver_impl(graph_impl_type& graph, const resolver_type& resolver)
                : base_type(graph, resolver){};

            const data_array_type& resolve1(const data_array_type& elevation)
            {
                return elevation;
            };
            const data_array_type& resolve2(const data_array_type& elevation)
            {
                return elevation;
            };
        };
    }

    struct no_sink_resolver
    {
        // TODO: flow graph implementation agnostic
        using flow_graph_impl_tag = detail::flow_graph_fixed_array_tag;
    };


    struct pflood_sink_resolver
    {
        // TODO: flow graph implementation agnostic
        using flow_graph_impl_tag = detail::flow_graph_fixed_array_tag;
    };


    namespace detail
    {

        template <class FG>
        class sink_resolver_impl<FG, pflood_sink_resolver>
            : public sink_resolver_impl_base<FG, pflood_sink_resolver>
        {
        public:
            using graph_impl_type = FG;
            using base_type = sink_resolver_impl_base<graph_impl_type, pflood_sink_resolver>;

            using data_array_type = typename graph_impl_type::data_array_type;

            sink_resolver_impl(graph_impl_type& graph, const pflood_sink_resolver& resolver)
                : base_type(graph, resolver){};

            const data_array_type& resolve1(const data_array_type& elevation)
            {
                m_filled_elevation = elevation;

                fill_sinks_sloped(this->m_graph_impl.grid(), m_filled_elevation);

                return m_filled_elevation;
            };

            const data_array_type& resolve2(const data_array_type& elevation)
            {
                return elevation;
            };

        private:
            data_array_type m_filled_elevation;
        };
    }


    enum class mst_route_method
    {
        basic,
        carve
    };


    struct mst_sink_resolver
    {
        using flow_graph_impl_tag = detail::flow_graph_fixed_array_tag;
        mst_method basin_method = mst_method::kruskal;
        mst_route_method route_method = mst_route_method::carve;
    };


    namespace detail
    {

        template <class FG>
        class sink_resolver_impl<FG, mst_sink_resolver>
            : public sink_resolver_impl_base<FG, mst_sink_resolver>
        {
        public:
            using graph_impl_type = FG;
            using base_type = sink_resolver_impl_base<graph_impl_type, mst_sink_resolver>;

            using size_type = typename graph_impl_type::size_type;
            using data_type = typename graph_impl_type::data_type;
            using data_array_type = typename graph_impl_type::data_array_type;

            sink_resolver_impl(graph_impl_type& graph, const mst_sink_resolver& resolver)
                : base_type(graph, resolver)
                , m_basin_graph(graph, resolver.basin_method){};

            const data_array_type& resolve1(const data_array_type& elevation)
            {
                return elevation;
            };

            const data_array_type& resolve2(const data_array_type& elevation)
            {
                // make sure the basins are up-to-date
                this->m_graph_impl.compute_basins();

                if (this->m_graph_impl.pits().empty())
                {
                    // return early, no sink to resolve
                    return elevation;
                }

                m_basin_graph.update_routes(elevation);

                if (this->m_resolver.route_method == mst_route_method::basic)
                {
                    update_routes_sinks_basic(elevation);
                }
                else
                {
                    // carve
                    update_routes_sinks_carve();
                }

                // finalize flow route update (donors and dfs graph traversal indices)
                this->m_graph_impl.compute_donors();
                this->m_graph_impl.compute_dfs_indices_downup();

                // fill sinks with tiny tilted surface
                fill_sinks_sloped(elevation);

                return m_filled_elevation;
            };

        private:
            basin_graph<FG> m_basin_graph;

            data_array_type m_filled_elevation;

            // basin graph edges are oriented in the counter flow direction
            static constexpr std::uint8_t outflow = 0;
            static constexpr std::uint8_t inflow = 1;

            void update_routes_sinks_basic(const data_array_type& elevation);
            void update_routes_sinks_carve();

            void fill_sinks_sloped(const data_array_type& elevation);
        };


        template <class FG>
        void sink_resolver_impl<FG, mst_sink_resolver>::update_routes_sinks_basic(
            const data_array_type& elevation)
        {
            auto& receivers = this->m_graph_impl.m_receivers;
            auto& dist2receivers = this->m_graph_impl.m_receivers_distance;

            // pits corresponds to outlets in inner basins
            const auto& pits = m_basin_graph.outlets();

            for (size_type edge_idx : m_basin_graph.tree())
            {
                auto& edge = m_basin_graph.edges()[edge_idx];

                // skip outer basins
                if (edge.pass[outflow] == static_cast<size_type>(-1))
                {
                    continue;
                }

                size_type pit_inflow = pits[edge.link[inflow]];

                // set infinite distance from the pit node to its receiver
                // (i.e., one of the two pass nodes).
                dist2receivers(pit_inflow, 0) = std::numeric_limits<data_type>::max();

                if (elevation.flat(edge.pass[inflow]) < elevation.flat(edge.pass[outflow]))
                {
                    receivers(pit_inflow, 0) = edge.pass[outflow];
                }
                else
                {
                    // also need to resolve the flow between the two
                    // nodes forming the pass, i.e., route the flow like this:
                    // pit -> pass (inflow) -> pass (outflow)
                    receivers(pit_inflow, 0) = edge.pass[inflow];
                    receivers(edge.pass[inflow], 0) = edge.pass[outflow];
                    dist2receivers(edge.pass[inflow], 0) = edge.pass_length;
                }
            }
        }


        template <class FG>
        void sink_resolver_impl<FG, mst_sink_resolver>::update_routes_sinks_carve()
        {
            auto& receivers = this->m_graph_impl.m_receivers;
            auto& dist2receivers = this->m_graph_impl.m_receivers_distance;

            // pits corresponds to outlets in inner basins
            const auto& pits = m_basin_graph.outlets();

            for (size_type edge_idx : m_basin_graph.tree())
            {
                auto& edge = m_basin_graph.edges()[edge_idx];

                // skip outer basins
                if (edge.pass[outflow] == static_cast<size_type>(-1))
                {
                    continue;
                }

                size_type pit_inflow = pits[edge.link[inflow]];

                // start at the pass (inflow) and then below follow the
                // receivers until the pit is found
                size_type current_node = edge.pass[inflow];
                size_type next_node = receivers(current_node, 0);
                data_type previous_dist = dist2receivers(current_node, 0);

                // re-route the pass inflow to the pass outflow
                receivers(current_node, 0) = edge.pass[outflow];
                dist2receivers(current_node, 0) = edge.pass_length;

                while (current_node != pit_inflow)
                {
                    auto rec_next_node = receivers(next_node, 0);
                    receivers(next_node, 0) = current_node;
                    std::swap(dist2receivers(next_node, 0), previous_dist);
                    current_node = next_node;
                    next_node = rec_next_node;
                }
            }
        }


        template <class FG>
        void sink_resolver_impl<FG, mst_sink_resolver>::fill_sinks_sloped(
            const data_array_type& elevation)
        {
            m_filled_elevation = elevation;

            const auto& dfs_indices = this->m_graph_impl.dfs_indices();
            const auto& receivers = this->m_graph_impl.receivers();

            for (const auto& idfs : dfs_indices)
            {
                const auto& irec = receivers(idfs, 0);

                if (idfs == irec)
                {
                    continue;
                }

                const auto& irec_elev = m_filled_elevation.flat(irec);

                if (m_filled_elevation.flat(idfs) <= m_filled_elevation.flat(irec))
                {
                    auto tiny_step
                        = std::nextafter(irec_elev, std::numeric_limits<data_type>::infinity());
                    m_filled_elevation.flat(idfs) = tiny_step;
                }
            }
        }
    }
}

#endif
